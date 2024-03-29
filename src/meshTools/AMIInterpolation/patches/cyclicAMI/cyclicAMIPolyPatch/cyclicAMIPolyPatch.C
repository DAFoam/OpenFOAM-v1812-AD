/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2018 OpenCFD Ltd.
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

#include "cyclicAMIPolyPatch.H"
#include "transformField.H"
#include "SubField.H"
#include "polyMesh.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"
#include "faceAreaIntersect.H"
#include "ops.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cyclicAMIPolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, cyclicAMIPolyPatch, word);
    addToRunTimeSelectionTable(polyPatch, cyclicAMIPolyPatch, dictionary);
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::vector Foam::cyclicAMIPolyPatch::findFaceNormalMaxRadius
(
    const pointField& faceCentres
) const
{
    // Determine a face furthest away from the axis

    const vectorField n((faceCentres - rotationCentre_) ^ rotationAxis_);

    const scalarField magRadSqr(magSqr(n));

    label facei = findMax(magRadSqr);

    if (debug)
    {
        Info<< "findFaceMaxRadius(const pointField&) : patch: " << name() << nl
            << "    rotFace  = " << facei << nl
            << "    point    = " << faceCentres[facei] << nl
            << "    distance = " << sqrt(magRadSqr[facei])
            << endl;
    }

    return n[facei];
}


void Foam::cyclicAMIPolyPatch::calcTransforms
(
    const primitivePatch& half0,
    const pointField& half0Ctrs,
    const vectorField& half0Areas,
    const pointField& half1Ctrs,
    const vectorField& half1Areas
)
{
    if (transform() != neighbPatch().transform())
    {
        FatalErrorInFunction
            << "Patch " << name()
            << " has transform type " << transformTypeNames[transform()]
            << ", neighbour patch " << neighbPatchName()
            << " has transform type "
            << neighbPatch().transformTypeNames[neighbPatch().transform()]
            << exit(FatalError);
    }


    // Calculate transformation tensors

    switch (transform())
    {
        case ROTATIONAL:
        {
            tensor revT = tensor::zero;

            if (rotationAngleDefined_)
            {
                const tensor T(rotationAxis_*rotationAxis_);

                const tensor S
                (
                    0, -rotationAxis_.z(), rotationAxis_.y(),
                    rotationAxis_.z(), 0, -rotationAxis_.x(),
                    -rotationAxis_.y(), rotationAxis_.x(), 0
                );

                const tensor revTPos
                (
                    T
                  + cos(rotationAngle_)*(tensor::I - T)
                  + sin(rotationAngle_)*S
                );

                const tensor revTNeg
                (
                    T
                  + cos(-rotationAngle_)*(tensor::I - T)
                  + sin(-rotationAngle_)*S
                );

                // Check - assume correct angle when difference in face areas
                // is the smallest
                const vector transformedAreaPos = gSum(half1Areas & revTPos);
                const vector transformedAreaNeg = gSum(half1Areas & revTNeg);
                const vector area0 = gSum(half0Areas);
                const scalar magArea0 = mag(area0) + ROOTVSMALL;

                // Areas have opposite sign, so sum should be zero when correct
                // rotation applied
                const scalar errorPos = mag(transformedAreaPos + area0);
                const scalar errorNeg = mag(transformedAreaNeg + area0);

                const scalar normErrorPos = errorPos/magArea0;
                const scalar normErrorNeg = errorNeg/magArea0;

                if (errorPos > errorNeg && normErrorNeg < matchTolerance())
                {
                    revT = revTNeg;
                    rotationAngle_ *= -1;
                }
                else
                {
                    revT = revTPos;
                }

                const scalar areaError = min(normErrorPos, normErrorNeg);

                if (areaError > matchTolerance())
                {
                    WarningInFunction
                        << "Patch areas are not consistent within "
                        << 100*matchTolerance()
                        << " % indicating a possible error in the specified "
                        << "angle of rotation" << nl
                        << "    owner patch     : " << name() << nl
                        << "    neighbour patch : " << neighbPatch().name()
                        << nl
                        << "    angle           : "
                        << radToDeg(rotationAngle_) << " deg" << nl
                        << "    area error      : " << 100*areaError << " %"
                        << "    match tolerance : " <<  matchTolerance()
                        << endl;
                }

                if (debug)
                {
                    scalar theta = radToDeg(rotationAngle_);

                    Pout<< "cyclicAMIPolyPatch::calcTransforms: patch:"
                        << name()
                        << " Specified rotation:"
                        << " swept angle: " << theta << " [deg]"
                        << " reverse transform: " << revT
                        << endl;
                }
            }
            else
            {
                point n0 = vector::zero;
                point n1 = vector::zero;
                if (half0Ctrs.size())
                {
                    n0 = findFaceNormalMaxRadius(half0Ctrs);
                }
                if (half1Ctrs.size())
                {
                    n1 = -findFaceNormalMaxRadius(half1Ctrs);
                }

                reduce(n0, maxMagSqrOp<point>());
                reduce(n1, maxMagSqrOp<point>());

                n0.normalise();
                n1.normalise();

                // Extended tensor from two local coordinate systems calculated
                // using normal and rotation axis
                const tensor E0
                (
                    rotationAxis_,
                    (n0 ^ rotationAxis_),
                    n0
                );
                const tensor E1
                (
                    rotationAxis_,
                    (-n1 ^ rotationAxis_),
                    -n1
                );
                revT = E1.T() & E0;

                if (debug)
                {
                    scalar theta = radToDeg(acos(-(n0 & n1)));

                    Pout<< "cyclicAMIPolyPatch::calcTransforms: patch:"
                        << name()
                        << " Specified rotation:"
                        << " n0:" << n0 << " n1:" << n1
                        << " swept angle: " << theta << " [deg]"
                        << " reverse transform: " << revT
                        << endl;
                }
            }

            const_cast<tensorField&>(forwardT()) = tensorField(1, revT.T());
            const_cast<tensorField&>(reverseT()) = tensorField(1, revT);
            const_cast<vectorField&>(separation()).setSize(0);
            const_cast<boolList&>(collocated()) = boolList(1, false);

            break;
        }
        case TRANSLATIONAL:
        {
            if (debug)
            {
                Pout<< "cyclicAMIPolyPatch::calcTransforms : patch:" << name()
                    << " Specified translation : " << separationVector_
                    << endl;
            }

            const_cast<tensorField&>(forwardT()).clear();
            const_cast<tensorField&>(reverseT()).clear();
            const_cast<vectorField&>(separation()) = vectorField
            (
                1,
                separationVector_
            );
            const_cast<boolList&>(collocated()) = boolList(1, false);

            break;
        }
        default:
        {
            if (debug)
            {
                Pout<< "patch:" << name()
                    << " Assuming cyclic AMI pairs are colocated" << endl;
            }

            const_cast<tensorField&>(forwardT()).clear();
            const_cast<tensorField&>(reverseT()).clear();
            const_cast<vectorField&>(separation()).setSize(0);
            const_cast<boolList&>(collocated()) = boolList(1, true);

            break;
        }
    }

    if (debug)
    {
        Pout<< "patch: " << name() << nl
            << "    forwardT = " << forwardT() << nl
            << "    reverseT = " << reverseT() << nl
            << "    separation = " << separation() << nl
            << "    collocated = " << collocated() << nl << endl;
    }
}


// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

void Foam::cyclicAMIPolyPatch::resetAMI
(
    const AMIPatchToPatchInterpolation::interpolationMethod& AMIMethod
) const
{
    if (owner())
    {
        AMIPtr_.clear();

        const polyPatch& nbr = neighbPatch();
        pointField nbrPoints
        (
            neighbPatch().boundaryMesh().mesh().points(),
            neighbPatch().meshPoints()
        );

        if (debug)
        {
            const Time& t = boundaryMesh().mesh().time();
            OFstream os(t.path()/name() + "_neighbourPatch-org.obj");
            meshTools::writeOBJ(os, neighbPatch().localFaces(), nbrPoints);
        }

        // Transform neighbour patch to local system
        transformPosition(nbrPoints);
        primitivePatch nbrPatch0
        (
            SubList<face>
            (
                nbr.localFaces(),
                nbr.size()
            ),
            nbrPoints
        );

        if (debug)
        {
            const Time& t = boundaryMesh().mesh().time();
            OFstream osN(t.path()/name() + "_neighbourPatch-trans.obj");
            meshTools::writeOBJ(osN, nbrPatch0.localFaces(), nbrPoints);

            OFstream osO(t.path()/name() + "_ownerPatch.obj");
            meshTools::writeOBJ(osO, this->localFaces(), localPoints());
        }

        // Construct/apply AMI interpolation to determine addressing and weights
        AMIPtr_.reset
        (
            new AMIPatchToPatchInterpolation
            (
                *this,
                nbrPatch0,
                surfPtr(),
                faceAreaIntersect::tmMesh,
                AMIRequireMatch_,
                AMIMethod,
                AMILowWeightCorrection_,
                AMIReverse_
            )
        );

        if (debug)
        {
            Pout<< "cyclicAMIPolyPatch : " << name()
                << " constructed AMI with " << nl
                << "    " << "srcAddress:" << AMIPtr_().srcAddress().size()
                << nl
                << "    " << "tgAddress :" << AMIPtr_().tgtAddress().size()
                << nl << endl;
        }
    }
}


void Foam::cyclicAMIPolyPatch::calcTransforms()
{
    const cyclicAMIPolyPatch& half0 = *this;
    vectorField half0Areas(half0.size());
    forAll(half0, facei)
    {
        half0Areas[facei] = half0[facei].areaNormal(half0.points());
    }

    const cyclicAMIPolyPatch& half1 = neighbPatch();
    vectorField half1Areas(half1.size());
    forAll(half1, facei)
    {
        half1Areas[facei] = half1[facei].areaNormal(half1.points());
    }

    calcTransforms
    (
        half0,
        half0.faceCentres(),
        half0Areas,
        half1.faceCentres(),
        half1Areas
    );

    if (debug)
    {
        Pout<< "calcTransforms() : patch: " << name() << nl
            << "    forwardT = " << forwardT() << nl
            << "    reverseT = " << reverseT() << nl
            << "    separation = " << separation() << nl
            << "    collocated = " << collocated() << nl << endl;
    }
}


void Foam::cyclicAMIPolyPatch::initGeometry(PstreamBuffers& pBufs)
{
    // The AMI is no longer valid. Leave it up to demand-driven calculation
    AMIPtr_.clear();

    polyPatch::initGeometry(pBufs);

    // Early calculation of transforms so e.g. cyclicACMI can use them.
    // Note: also triggers primitiveMesh face centre. Note that cell
    // centres should -not- be calculated
    // since e.g. cyclicACMI override face areas
    calcTransforms();
}


void Foam::cyclicAMIPolyPatch::calcGeometry(PstreamBuffers& pBufs)
{
    // All geometry done inside initGeometry
}


void Foam::cyclicAMIPolyPatch::initMovePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    // The AMI is no longer valid. Leave it up to demand-driven calculation
    AMIPtr_.clear();

    polyPatch::initMovePoints(pBufs, p);

    // See below. Clear out any local geometry
    primitivePatch::movePoints(p);

    // Early calculation of transforms. See above.
    calcTransforms();
}


void Foam::cyclicAMIPolyPatch::movePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    polyPatch::movePoints(pBufs, p);

    // All transformation tensors already done in initMovePoints
}


void Foam::cyclicAMIPolyPatch::initUpdateMesh(PstreamBuffers& pBufs)
{
    // The AMI is no longer valid. Leave it up to demand-driven calculation
    AMIPtr_.clear();

    polyPatch::initUpdateMesh(pBufs);
}


void Foam::cyclicAMIPolyPatch::updateMesh(PstreamBuffers& pBufs)
{
    polyPatch::updateMesh(pBufs);
}


void Foam::cyclicAMIPolyPatch::clearGeom()
{
    AMIPtr_.clear();
    polyPatch::clearGeom();
}


// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * //

Foam::cyclicAMIPolyPatch::cyclicAMIPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType,
    const transformType transform
)
:
    coupledPolyPatch(name, size, start, index, bm, patchType, transform),
    nbrPatchName_(word::null),
    nbrPatchID_(-1),
    rotationAxis_(vector::zero),
    rotationCentre_(vector::zero),
    rotationAngleDefined_(false),
    rotationAngle_(0.0),
    separationVector_(vector::zero),
    AMIPtr_(nullptr),
    AMIMethod_(AMIPatchToPatchInterpolation::imFaceAreaWeight),
    AMIReverse_(false),
    AMIRequireMatch_(true),
    AMILowWeightCorrection_(-1.0),
    surfPtr_(nullptr),
    surfDict_(fileName("surface"))
{
    // Neighbour patch might not be valid yet so no transformation
    // calculation possible
}


Foam::cyclicAMIPolyPatch::cyclicAMIPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    coupledPolyPatch(name, dict, index, bm, patchType),
    nbrPatchName_(dict.lookupOrDefault<word>("neighbourPatch", "")),
    coupleGroup_(dict),
    nbrPatchID_(-1),
    rotationAxis_(vector::zero),
    rotationCentre_(vector::zero),
    rotationAngleDefined_(false),
    rotationAngle_(0.0),
    separationVector_(vector::zero),
    AMIPtr_(nullptr),
    AMIMethod_
    (
        AMIPatchToPatchInterpolation::interpolationMethodNames_
        [
            dict.lookupOrDefault
            (
                "method",
                AMIPatchToPatchInterpolation::interpolationMethodNames_
                [
                    AMIPatchToPatchInterpolation::imFaceAreaWeight
                ]
            )
        ]
    ),
    AMIReverse_(dict.lookupOrDefault("flipNormals", false)),
    AMIRequireMatch_(true),
    AMILowWeightCorrection_(dict.lookupOrDefault("lowWeightCorrection", -1.0)),
    surfPtr_(nullptr),
    surfDict_(dict.subOrEmptyDict("surface"))
{
    if (nbrPatchName_ == word::null && !coupleGroup_.valid())
    {
        FatalIOErrorInFunction(dict)
            << "No \"neighbourPatch\" or \"coupleGroup\" provided."
            << exit(FatalIOError);
    }

    if (nbrPatchName_ == name)
    {
        FatalIOErrorInFunction(dict)
            << "Neighbour patch name " << nbrPatchName_
            << " cannot be the same as this patch " << name
            << exit(FatalIOError);
    }

    switch (transform())
    {
        case ROTATIONAL:
        {
            dict.readEntry("rotationAxis", rotationAxis_);
            dict.readEntry("rotationCentre", rotationCentre_);
            if (dict.readIfPresent("rotationAngle", rotationAngle_))
            {
                rotationAngleDefined_ = true;
                rotationAngle_ = degToRad(rotationAngle_);

                if (debug)
                {
                    Info<< "rotationAngle: " << rotationAngle_ << " [rad]"
                        <<  endl;
                }
            }

            scalar magRot = mag(rotationAxis_);
            if (magRot < SMALL)
            {
                FatalIOErrorInFunction(dict)
                    << "Illegal rotationAxis " << rotationAxis_ << endl
                    << "Please supply a non-zero vector."
                    << exit(FatalIOError);
            }
            rotationAxis_ /= magRot;

            break;
        }
        case TRANSLATIONAL:
        {
            dict.readEntry("separationVector", separationVector_);
            break;
        }
        default:
        {
            // No additional info required
        }
    }

    // Neighbour patch might not be valid yet so no transformation
    // calculation possible
}


Foam::cyclicAMIPolyPatch::cyclicAMIPolyPatch
(
    const cyclicAMIPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    coupledPolyPatch(pp, bm),
    nbrPatchName_(pp.nbrPatchName_),
    coupleGroup_(pp.coupleGroup_),
    nbrPatchID_(-1),
    rotationAxis_(pp.rotationAxis_),
    rotationCentre_(pp.rotationCentre_),
    rotationAngleDefined_(pp.rotationAngleDefined_),
    rotationAngle_(pp.rotationAngle_),
    separationVector_(pp.separationVector_),
    AMIPtr_(nullptr),
    AMIMethod_(pp.AMIMethod_),
    AMIReverse_(pp.AMIReverse_),
    AMIRequireMatch_(pp.AMIRequireMatch_),
    AMILowWeightCorrection_(pp.AMILowWeightCorrection_),
    surfPtr_(nullptr),
    surfDict_(pp.surfDict_)
{
    // Neighbour patch might not be valid yet so no transformation
    // calculation possible
}


Foam::cyclicAMIPolyPatch::cyclicAMIPolyPatch
(
    const cyclicAMIPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart,
    const word& nbrPatchName
)
:
    coupledPolyPatch(pp, bm, index, newSize, newStart),
    nbrPatchName_(nbrPatchName),
    coupleGroup_(pp.coupleGroup_),
    nbrPatchID_(-1),
    rotationAxis_(pp.rotationAxis_),
    rotationCentre_(pp.rotationCentre_),
    rotationAngleDefined_(pp.rotationAngleDefined_),
    rotationAngle_(pp.rotationAngle_),
    separationVector_(pp.separationVector_),
    AMIPtr_(nullptr),
    AMIMethod_(pp.AMIMethod_),
    AMIReverse_(pp.AMIReverse_),
    AMIRequireMatch_(pp.AMIRequireMatch_),
    AMILowWeightCorrection_(pp.AMILowWeightCorrection_),
    surfPtr_(nullptr),
    surfDict_(pp.surfDict_)
{
    if (nbrPatchName_ == name())
    {
        FatalErrorInFunction
            << "Neighbour patch name " << nbrPatchName_
            << " cannot be the same as this patch " << name()
            << exit(FatalError);
    }

    // Neighbour patch might not be valid yet so no transformation
    // calculation possible
}


Foam::cyclicAMIPolyPatch::cyclicAMIPolyPatch
(
    const cyclicAMIPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const label newStart
)
:
    coupledPolyPatch(pp, bm, index, mapAddressing, newStart),
    nbrPatchName_(pp.nbrPatchName_),
    coupleGroup_(pp.coupleGroup_),
    nbrPatchID_(-1),
    rotationAxis_(pp.rotationAxis_),
    rotationCentre_(pp.rotationCentre_),
    rotationAngleDefined_(pp.rotationAngleDefined_),
    rotationAngle_(pp.rotationAngle_),
    separationVector_(pp.separationVector_),
    AMIPtr_(nullptr),
    AMIMethod_(pp.AMIMethod_),
    AMIReverse_(pp.AMIReverse_),
    AMIRequireMatch_(pp.AMIRequireMatch_),
    AMILowWeightCorrection_(pp.AMILowWeightCorrection_),
    surfPtr_(nullptr),
    surfDict_(pp.surfDict_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cyclicAMIPolyPatch::~cyclicAMIPolyPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::cyclicAMIPolyPatch::neighbPatchID() const
{
    if (nbrPatchID_ == -1)
    {
        nbrPatchID_ = this->boundaryMesh().findPatchID(neighbPatchName());

        if (nbrPatchID_ == -1)
        {
            FatalErrorInFunction
                << "Illegal neighbourPatch name " << neighbPatchName()
                << nl << "Valid patch names are "
                << this->boundaryMesh().names()
                << exit(FatalError);
        }

        // Check that it is a cyclic AMI patch
        const cyclicAMIPolyPatch& nbrPatch =
            refCast<const cyclicAMIPolyPatch>
            (
                this->boundaryMesh()[nbrPatchID_]
            );

        if (nbrPatch.neighbPatchName() != name())
        {
            WarningInFunction
                << "Patch " << name()
                << " specifies neighbour patch " << neighbPatchName()
                << nl << " but that in return specifies "
                << nbrPatch.neighbPatchName() << endl;
        }
    }

    return nbrPatchID_;
}


bool Foam::cyclicAMIPolyPatch::owner() const
{
    return index() < neighbPatchID();
}


const Foam::cyclicAMIPolyPatch& Foam::cyclicAMIPolyPatch::neighbPatch() const
{
    const polyPatch& pp = this->boundaryMesh()[neighbPatchID()];
    return refCast<const cyclicAMIPolyPatch>(pp);
}


const Foam::autoPtr<Foam::searchableSurface>&
Foam::cyclicAMIPolyPatch::surfPtr() const
{
    const word surfType(surfDict_.lookupOrDefault<word>("type", "none"));

    if (!surfPtr_.valid() && owner() && surfType != "none")
    {
        word surfName(surfDict_.lookupOrDefault("name", name()));

        const polyMesh& mesh = boundaryMesh().mesh();

        surfPtr_ =
            searchableSurface::New
            (
                surfType,
                IOobject
                (
                    surfName,
                    mesh.time().constant(),
                    "triSurface",
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                surfDict_
            );
    }

    return surfPtr_;
}


const Foam::AMIPatchToPatchInterpolation& Foam::cyclicAMIPolyPatch::AMI() const
{
    if (!owner())
    {
        FatalErrorInFunction
            << "AMI interpolator only available to owner patch"
            << abort(FatalError);
    }

    if (!AMIPtr_.valid())
    {
        resetAMI(AMIMethod_);
    }

    return *AMIPtr_;
}


bool Foam::cyclicAMIPolyPatch::applyLowWeightCorrection() const
{
    if (owner())
    {
        return AMI().applyLowWeightCorrection();
    }
    else
    {
        return neighbPatch().AMI().applyLowWeightCorrection();
    }
}


void Foam::cyclicAMIPolyPatch::transformPosition(pointField& l) const
{
    if (!parallel())
    {
        if (transform() == ROTATIONAL)
        {
            l = Foam::transform(forwardT(), l - rotationCentre_)
              + rotationCentre_;
        }
        else
        {
            l = Foam::transform(forwardT(), l);
        }
    }
    else if (separated())
    {
        // transformPosition gets called on the receiving side,
        // separation gets calculated on the sending side so subtract

        const vectorField& s = separation();
        if (s.size() == 1)
        {
            forAll(l, i)
            {
                l[i] -= s[0];
            }
        }
        else
        {
            l -= s;
        }
    }
}


void Foam::cyclicAMIPolyPatch::transformPosition
(
    point& l,
    const label facei
) const
{
    if (!parallel())
    {
        const tensor& T =
        (
            forwardT().size() == 1
          ? forwardT()[0]
          : forwardT()[facei]
        );

        if (transform() == ROTATIONAL)
        {
            l = Foam::transform(T, l - rotationCentre_) + rotationCentre_;
        }
        else
        {
            l = Foam::transform(T, l);
        }
    }
    else if (separated())
    {
        const vector& s =
        (
            separation().size() == 1
          ? separation()[0]
          : separation()[facei]
        );

        l -= s;
    }
}


void Foam::cyclicAMIPolyPatch::reverseTransformPosition
(
    point& l,
    const label facei
) const
{
    if (!parallel())
    {
        const tensor& T =
        (
            reverseT().size() == 1
          ? reverseT()[0]
          : reverseT()[facei]
        );

        if (transform() == ROTATIONAL)
        {
            l = Foam::transform(T, l - rotationCentre_) + rotationCentre_;
        }
        else
        {
            l = Foam::transform(T, l);
        }
    }
    else if (separated())
    {
        const vector& s =
        (
            separation().size() == 1
          ? separation()[0]
          : separation()[facei]
        );

        l += s;
    }
}


void Foam::cyclicAMIPolyPatch::reverseTransformDirection
(
    vector& d,
    const label facei
) const
{
    if (!parallel())
    {
        const tensor& T =
        (
            reverseT().size() == 1
          ? reverseT()[0]
          : reverseT()[facei]
        );

        d = Foam::transform(T, d);
    }
}


void Foam::cyclicAMIPolyPatch::calcGeometry
(
    const primitivePatch& referPatch,
    const pointField& thisCtrs,
    const vectorField& thisAreas,
    const pointField& thisCc,
    const pointField& nbrCtrs,
    const vectorField& nbrAreas,
    const pointField& nbrCc
)
{}


void Foam::cyclicAMIPolyPatch::initOrder
(
    PstreamBuffers& pBufs,
    const primitivePatch& pp
) const
{}


bool Foam::cyclicAMIPolyPatch::order
(
    PstreamBuffers& pBufs,
    const primitivePatch& pp,
    labelList& faceMap,
    labelList& rotation
) const
{
    faceMap.setSize(pp.size());
    faceMap = -1;

    rotation.setSize(pp.size());
    rotation = 0;

    return false;
}


Foam::label Foam::cyclicAMIPolyPatch::pointFace
(
    const label facei,
    const vector& n,
    point& p
) const
{
    point prt(p);
    reverseTransformPosition(prt, facei);

    vector nrt(n);
    reverseTransformDirection(nrt, facei);

    label nbrFacei = -1;

    if (owner())
    {
        nbrFacei = AMI().tgtPointFace
        (
            *this,
            neighbPatch(),
            nrt,
            facei,
            prt
        );
    }
    else
    {
        nbrFacei = neighbPatch().AMI().srcPointFace
        (
            neighbPatch(),
            *this,
            nrt,
            facei,
            prt
        );
    }

    if (nbrFacei >= 0)
    {
        p = prt;
    }

    return nbrFacei;
}


void Foam::cyclicAMIPolyPatch::write(Ostream& os) const
{
    coupledPolyPatch::write(os);
    if (!nbrPatchName_.empty())
    {
        os.writeEntry("neighbourPatch", nbrPatchName_);
    }
    coupleGroup_.write(os);

    switch (transform())
    {
        case ROTATIONAL:
        {
            os.writeEntry("rotationAxis", rotationAxis_);
            os.writeEntry("rotationCentre", rotationCentre_);

            if (rotationAngleDefined_)
            {
                os.writeEntry("rotationAngle", radToDeg(rotationAngle_));
            }

            break;
        }
        case TRANSLATIONAL:
        {
            os.writeEntry("separationVector", separationVector_);
            break;
        }
        case NOORDERING:
        {
            break;
        }
        default:
        {
            // No additional info to write
        }
    }

    if (AMIMethod_ != AMIPatchToPatchInterpolation::imFaceAreaWeight)
    {
        os.writeEntry
        (
            "method",
            AMIPatchToPatchInterpolation::interpolationMethodNames_
            [
                AMIMethod_
            ]
        );
    }

    if (AMIReverse_)
    {
        os.writeEntry("flipNormals", AMIReverse_);
    }

    if (AMILowWeightCorrection_ > 0)
    {
        os.writeEntry("lowWeightCorrection", AMILowWeightCorrection_);
    }

    if (!surfDict_.empty())
    {
        surfDict_.writeEntry(surfDict_.dictName(), os);
    }
}


// ************************************************************************* //
