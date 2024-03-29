/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Author: Franjo Juretic (franjo.juretic@c-fields.com)
     \\/     M anipulation  | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "meshOctree.H"
#include "triSurf.H"
#include "demandDrivenData.H"
#include "helperFunctions.H"
#include "HashSet.H"

// #define DEBUGSearch

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::meshOctree::findNearestSurfacePoint
(
    point& nearest,
    scalar& distSq,
    label& nearestTriangle,
    label& region,
    const point& p
) const
{
    region = -1;
    nearestTriangle = 1;

    const label cLabel = findLeafContainingVertex(p);
    vector sizeVec;
    if (cLabel < 0)
    {
        sizeVec.x() = sizeVec.y() = sizeVec.z() = searchRange_;
    }
    else
    {
        const scalar s = 0.75*leaves_[cLabel]->size(rootBox_);
        sizeVec.x() = sizeVec.y() = sizeVec.z() = s;
    }

    // find nearest surface vertex to the point p
    bool found(false);
    label iterationI(0);
    DynList<const meshOctreeCube*, 256> neighbours;

    distSq = VGREAT;

    do
    {
        const boundBox bb(p - sizeVec, p + sizeVec);

        neighbours.clear();
        findLeavesContainedInBox(bb, neighbours);

        // find nearest projection
        forAll(neighbours, neiI)
        {
            if (!neighbours[neiI]->hasContainedElements())
                continue;

            const VRWGraph& ct =
                neighbours[neiI]->slotPtr()->containedTriangles_;
            const constRow el = ct[neighbours[neiI]->containedElements()];
            forAll(el, tI)
            {
                const point p0 =
                    help::nearestPointOnTheTriangle(el[tI], surface_, p);

                const scalar dSq = Foam::magSqr(p0 - p);
                if (dSq < distSq)
                {
                    distSq = dSq;
                    nearest = p0;
                    nearestTriangle = el[tI];
                    region = surface_[el[tI]].region();
                    found = true;
                }
            }
        }

        if (!found)
            sizeVec *= 2.0;

    } while (!found && (iterationI++< 100));

    # ifdef DEBUGSearch
    forAll(surface_, triI)
    {
        const point pp = help::nearestPointOnTheTriangle(triI, surface_, p);

        if (distSq - magSqr(pp - p) > SMALL)
            Pout<< "Point " << p << " current nearest " << nearest
                << " closer point " << pp << endl;
    }
    # endif

    if ((!found || (region < 0)) && !Pstream::parRun())
    {
        Warning << "Could not find a boundary region for vertex " << p << endl;
        Warning << "Found " << found << " and region " << region << endl;
    }
}


void Foam::Module::meshOctree::findNearestSurfacePointInRegion
(
    point& nearest,
    scalar& distSq,
    label& nearestTriangle,
    const label region,
    const point& p
) const
{
    if (region < 0)
    {
        WarningInFunction
            << "Region " << region << " is not valid!" << endl;

        return;
    }
    const label cLabel = findLeafContainingVertex(p);
    vector sizeVec;
    if (cLabel < 0)
    {
        sizeVec.x() = sizeVec.y() = sizeVec.z() = searchRange_;
    }
    else
    {
        const scalar s = 0.75*leaves_[cLabel]->size(rootBox_);
        sizeVec.x() = sizeVec.y() = sizeVec.z() = s;
    }

    // find nearest surface vertex to the point p
    bool found(false);
    label iterationI(0);
    DynList<const meshOctreeCube*, 256> neighbours;
    nearestTriangle = -1;

    distSq = VGREAT;

    do
    {
        const boundBox bb(p - sizeVec, p + sizeVec);

        neighbours.clear();
        findLeavesContainedInBox(bb, neighbours);

        // find nearest projection
        forAll(neighbours, neiI)
        {
            if (!neighbours[neiI]->hasContainedElements())
                continue;

            const VRWGraph& ct =
                neighbours[neiI]->slotPtr()->containedTriangles_;
            const constRow el =
                ct[neighbours[neiI]->containedElements()];
            forAll(el, tI)
            {
                if (surface_[el[tI]].region() != region)
                    continue;

                const point p0 =
                    help::nearestPointOnTheTriangle(el[tI], surface_, p);

                const scalar dSq = Foam::magSqr(p0 - p);
                if (dSq < distSq)
                {
                    distSq = dSq;
                    nearest = p0;
                    nearestTriangle = el[tI];
                    found = true;
                }
            }
        }

        if (!found)
            sizeVec *= 2.0;

    } while (!found && (iterationI++< 5));

    if ((!found || (region < 0)) && !Pstream::parRun())
        Warning << "Could not find a boundary region for vertex " << p << endl;
}


bool Foam::Module::meshOctree::findNearestEdgePoint
(
    point& edgePoint,
    scalar& distSq,
    label& nearestEdge,
    const point& p,
    const DynList<label>& regions
) const
{
    // find the estimate for the searching range
    const label cLabel = findLeafContainingVertex(p);
    vector sizeVec;
    if (cLabel < 0)
    {
        sizeVec.x() = sizeVec.y() = sizeVec.z() = searchRange_;
    }
    else
    {
        const scalar s = 0.75*leaves_[cLabel]->size(rootBox_);
        sizeVec.x() = sizeVec.y() = sizeVec.z() = s;
    }

    DynList<const meshOctreeCube*, 256> neighbours;

    const pointField& sPoints = surface_.points();
    const edgeLongList& edges = surface_.edges();
    const VRWGraph& edgeFaces = surface_.edgeFacets();

    edgePoint = p;
    bool foundAnEdge(false);
    label iterationI(0);

    distSq = VGREAT;
    nearestEdge = -1;

//    Info<< "Finding nearest point for " << p << " size vec " << sizeVec << endl;

    do
    {
        const boundBox bb(p - sizeVec, p + sizeVec);

        neighbours.clear();
        findLeavesContainedInBox(bb, neighbours);

//        Info<< "Iteration " << iterationI << " nu found boxes "
//             << neighbours.size() << endl;

        forAll(neighbours, neiI)
        {
            if (!neighbours[neiI]->hasContainedEdges())
                continue;

            const VRWGraph& containedEdges =
                neighbours[neiI]->slotPtr()->containedEdges_;
            const constRow ce =
                containedEdges[neighbours[neiI]->containedEdges()];

//            Info<< "Number of contained edges in box "
//                 << neighbours[neiI]->cubeLabel()
//                 << " are " << ce.size() << endl;

            forAll(ce, eI)
            {
                const label edgeI = ce[eI];

                // find if the edge is in correct patches
                bool correctPatches(true);

                forAllRow(edgeFaces, edgeI, efI)
                {
                    const label facetI = edgeFaces(edgeI, efI);

                    if (!regions.found(surface_[facetI].region()))
                    {
                        correctPatches = false;
                        break;
                    }
                }

                if (!correctPatches)
                    continue;

                const edge& e = edges[edgeI];
                const point& sp = sPoints[e.start()];
                const point& ep = sPoints[e.end()];
                const point np = help::nearestPointOnTheEdgeExact(sp, ep, p);
                const scalar dSq = Foam::magSqr(np - p);

                if (dSq < distSq)
                {
                    distSq = dSq;
                    edgePoint = np;
                    nearestEdge = ce[eI];
                    foundAnEdge = true;
                }
            }
        }

        if (!foundAnEdge)
            sizeVec *= 2.0;

    } while (!foundAnEdge && (++iterationI < 3));

    return foundAnEdge;
}


bool Foam::Module::meshOctree::findNearestPointToEdge
(
    point& nearest,
    scalar& distSq,
    label& nearestEdge,
    const FixedList<point, 2>& edgePoints,
    const FixedList<label, 2>& edgePointRegions
) const
{
    const point c = 0.5*(edgePoints[0] + edgePoints[1]);
    const scalar dst = mag(edgePoints[0] - edgePoints[1]);
    vector sizeVec(dst, dst, dst);

    boundBox bb(c - 0.75*sizeVec, c + 0.75*sizeVec);

    DynList<const meshOctreeCube*, 256> leavesInBox;
    findLeavesContainedInBox(bb, leavesInBox);

    const VRWGraph& edgeFaces = surface_.edgeFacets();
    const pointField& points = surface_.points();
    const edgeLongList& surfaceEdges = surface_.edges();

    distSq = VGREAT;
    nearestEdge = -1;

    bool found(false);

    forAll(leavesInBox, leafI)
    {
        if (!leavesInBox[leafI]->hasContainedEdges())
            continue;

        const VRWGraph& containedEdges =
            leavesInBox[leafI]->slotPtr()->containedEdges_;
        const constRow edges =
            containedEdges[leavesInBox[leafI]->containedEdges()];

        forAll(edges, eI)
        {
            const constRow ef = edgeFaces[edges[eI]];
            if (ef.size() != 2)
                continue;

            if
            (
                (
                    (surface_[ef[0]].region() == edgePointRegions[0])
                 && (surface_[ef[1]].region() == edgePointRegions[1])
                )
             || (
                    (surface_[ef[1]].region() == edgePointRegions[0])
                 && (surface_[ef[0]].region() == edgePointRegions[1])
                )
            )
            {
                const edge& edg = surfaceEdges[edges[eI]];

                point nearestOnEdge, nearestOnLine;
                if (
                    help::nearestEdgePointToTheLine
                    (
                        points[edg[0]],
                        points[edg[1]],
                        edgePoints[0],
                        edgePoints[1],
                        nearestOnEdge,
                        nearestOnLine
                    )
                )
                {
                    if (magSqr(nearestOnEdge-nearestOnLine) < distSq)
                    {
                        nearest = nearestOnEdge;
                        nearestEdge = edges[eI];
                        distSq = magSqr(nearestOnEdge-nearestOnLine);
                        found = true;
                    }
                }
            }
        }
    }

    return found;
}


bool Foam::Module::meshOctree::findNearestCorner
(
    point& nearest,
    scalar& distSq,
    label& nearestPoint,
    const point& p,
    const DynList<label>& patches
) const
{
    const label cLabel = findLeafContainingVertex(p);
    vector sizeVec;
    if (cLabel < 0)
    {
        sizeVec.x() = sizeVec.y() = sizeVec.z() = searchRange_;
    }
    else
    {
        const scalar s = 0.75*leaves_[cLabel]->size(rootBox_);
        sizeVec.x() = sizeVec.y() = sizeVec.z() = s;
    }

    // find nearest surface vertex to the point p
    bool found(false);
    label iterationI(0);
    DynList<const meshOctreeCube*, 256> neighbours;

    const pointField& points = surface_.points();
    const VRWGraph& pEdges = surface_.pointEdges();
    const VRWGraph& eFacets = surface_.edgeFacets();

    distSq = VGREAT;
    nearestPoint = -1;

    do
    {
        const boundBox bb(p - sizeVec, p + sizeVec);

        neighbours.clear();
        findLeavesContainedInBox(bb, neighbours);
        labelHashSet checkedPoint;

        // find nearest projection
        forAll(neighbours, neiI)
        {
            if (!neighbours[neiI]->hasContainedElements())
                continue;

            const VRWGraph& ct =
                neighbours[neiI]->slotPtr()->containedTriangles_;
            const constRow el = ct[neighbours[neiI]->containedElements()];
            forAll(el, tI)
            {
                const labelledTri& tri = surface_[el[tI]];

                forAll(tri, pI)
                {
                    const label spI = tri[pI];

                    if (checkedPoint.found(spI))
                        continue;

                    checkedPoint.insert(spI);

                    DynList<label> nodePatches;
                    label nEdges(0);

                    forAllRow(pEdges, spI, i)
                    {
                        const label eI = pEdges(spI, i);

                        if (eFacets.sizeOfRow(eI) != 2)
                            break;

                        if
                        (
                            surface_[eFacets(eI, 0)].region()
                         != surface_[eFacets(eI, 1)].region()
                        )
                        {
                            // found an edge attached to this vertex
                            ++nEdges;
                            nodePatches.appendIfNotIn
                            (
                                surface_[eFacets(eI, 0)].region()
                            );
                            nodePatches.appendIfNotIn
                            (
                                surface_[eFacets(eI, 1)].region()
                            );
                        }
                    }

                    if (nEdges > 2)
                    {
                        // check if all required patches
                        // are present at this corner
                        nEdges = 0;
                        forAll(patches, i)
                        {
                            if (nodePatches.found(patches[i]))
                            {
                                ++nEdges;
                            }
                        }

                        if (nEdges >= patches.size())
                        {
                            // all patches are present, check the distance
                            const scalar dSq = Foam::magSqr(points[spI] - p);

                            if (dSq < distSq)
                            {
                                distSq = dSq;
                                found = true;
                                nearest = points[spI];
                                nearestPoint = spI;
                            }
                        }
                    }
                }
            }
        }

        if (!found)
            sizeVec *= 2.0;

    } while (!found && (iterationI++< 3));

    return found;
}


bool Foam::Module::meshOctree::findNearestPointToPatches
(
    point& nearest,
    scalar& distSq,
    const point& p,
    const DynList<label>& patches,
    const scalar tol
) const
{
    if (patches.size() == 0)
        return false;

    nearest = p;
    scalar distSqApprox;
    label iter(0);
    while (iter++ < 40)
    {
        point newP(vector::zero);
        forAll(patches, patchI)
        {
            point np;
            label nearestTri;
            this->findNearestSurfacePointInRegion
            (
                np,
                distSqApprox,
                nearestTri,
                patches[patchI],
                nearest
            );

            newP += np;
        }

        newP /= patches.size();
        distSq = magSqr(newP - p);

        if (Foam::magSqr(newP - nearest) < tol*distSq)
            break;

        nearest = newP;
    }

    return true;
}


// ************************************************************************* //
