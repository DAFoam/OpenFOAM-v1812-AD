/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017-2018 OpenCFD Ltd.
     \\/     M anipulation  |
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

#include "vtkWrite.H"
#include "dictionary.H"
#include "Time.H"
#include "areaFields.H"
#include "stringListOps.H"
#include "foamVtkInternalWriter.H"
#include "foamVtkPatchWriter.H"
#include "foamVtkSeriesWriter.H"
#include "foamVtmWriter.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(vtkWrite, 0);
    addToRunTimeSelectionTable(functionObject, vtkWrite, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::functionObjects::vtkWrite::writeAllVolFields
(
    autoPtr<vtk::internalWriter>& internalWriter,
    UPtrList<vtk::patchWriter>& patchWriters,
    const fvMeshSubset& proxy,
    const wordHashSet& acceptField
) const
{
    #undef  vtkWrite_WRITE_FIELD
    #define vtkWrite_WRITE_FIELD(FieldType)     \
        writeVolFields<FieldType>               \
        (                                       \
            internalWriter,                     \
            patchWriters,                       \
            proxy,                              \
            acceptField                         \
        )


    label count = 0;
    count += vtkWrite_WRITE_FIELD(volScalarField);
    count += vtkWrite_WRITE_FIELD(volVectorField);
    count += vtkWrite_WRITE_FIELD(volSphericalTensorField);
    count += vtkWrite_WRITE_FIELD(volSymmTensorField);
    count += vtkWrite_WRITE_FIELD(volTensorField);

    #undef vtkWrite_WRITE_FIELD
    return count;
}


Foam::label Foam::functionObjects::vtkWrite::writeAllVolFields
(
    autoPtr<vtk::internalWriter>& internalWriter,
    const autoPtr<volPointInterpolation>& pInterp,

    UPtrList<vtk::patchWriter>& patchWriters,
    const UPtrList<PrimitivePatchInterpolation<primitivePatch>>& patchInterps,
    const fvMeshSubset& proxy,
    const wordHashSet& acceptField
) const
{
    #undef  vtkWrite_WRITE_FIELD
    #define vtkWrite_WRITE_FIELD(FieldType)     \
        writeVolFields<FieldType>               \
        (                                       \
            internalWriter, pInterp,            \
            patchWriters,   patchInterps,       \
            proxy,                              \
            acceptField                         \
        )


    label count = 0;
    count += vtkWrite_WRITE_FIELD(volScalarField);
    count += vtkWrite_WRITE_FIELD(volVectorField);
    count += vtkWrite_WRITE_FIELD(volSphericalTensorField);
    count += vtkWrite_WRITE_FIELD(volSymmTensorField);
    count += vtkWrite_WRITE_FIELD(volTensorField);

    #undef vtkWrite_WRITE_FIELD
    return count;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::vtkWrite::vtkWrite
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObject(name),
    time_(runTime),
    outputDir_(),
    printf_(),
    writeOpts_(vtk::formatType::INLINE_BASE64),
    verbose_(true),
    doInternal_(true),
    doBoundary_(true),
    oneBoundary_(false),
    interpolate_(false),
    decompose_(false),
    writeIds_(false),
    meshState_(polyMesh::TOPO_CHANGE),
    selectRegions_(),
    selectPatches_(),
    selectFields_(),
    selection_(),
    meshes_(),
    meshSubsets_(),
    vtuMappings_(),
    series_()
{
    // May still want this? (OCT-2018)
    // if (postProcess)
    // {
    //     // Disable for post-process mode.
    //     // Emit as FatalError for the try/catch in the caller.
    //     FatalError
    //         << type() << " disabled in post-process mode"
    //         << exit(FatalError);
    // }

    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::vtkWrite::read(const dictionary& dict)
{
    functionObject::read(dict);

    readSelection(dict);

    // We probably cannot trust old information after a reread
    series_.clear();

    // verbose_ = dict.lookupOrDefault<bool>("verbose", true);
    doInternal_ = dict.lookupOrDefault<bool>("internal", true);
    doBoundary_ = dict.lookupOrDefault<bool>("boundary", true);
    oneBoundary_ = dict.lookupOrDefault<bool>("single", false);
    interpolate_ = dict.lookupOrDefault<bool>("interpolate", false);

    //
    // Writer options - default is xml base64
    //
    writeOpts_ = vtk::formatType::INLINE_BASE64;

    writeOpts_.ascii
    (
        dict.found("format")
     && (IOstream::formatEnum(dict.get<word>("format")) == IOstream::ASCII)
    );

    if (dict.lookupOrDefault("legacy", false))
    {
        writeOpts_.legacy(true);
    }

    writeOpts_.precision
    (
        dict.lookupOrDefault
        (
            "precision",
            IOstream::defaultPrecision()
        )
    );

    // Info<< type() << " " << name() << " output-format: "
    //     << writeOpts_.description() << nl;

    const int padWidth = dict.lookupOrDefault<int>("width", 8);

    // Appropriate printf format - Enforce min/max sanity limits
    if (padWidth < 1 || padWidth > 31)
    {
        printf_.clear();
    }
    else
    {
        printf_ = "%0" + std::to_string(padWidth) + "d";
    }

    //
    // Other options
    //

    decompose_ = dict.lookupOrDefault("decompose", false);
    writeIds_ = dict.lookupOrDefault("writeIds", false);


    // Output directory

    outputDir_.clear();
    dict.readIfPresent("directory", outputDir_);

    if (outputDir_.size())
    {
        // User-defined output directory
        outputDir_.expand();
        if (!outputDir_.isAbsolute())
        {
            outputDir_ = time_.globalPath()/outputDir_;
        }
    }
    else
    {
        // Standard postProcessing/ naming
        outputDir_ = time_.globalPath()/functionObject::outputPrefix/name();
    }
    outputDir_.clean();

    return true;
}


bool Foam::functionObjects::vtkWrite::execute()
{
    return true;
}


bool Foam::functionObjects::vtkWrite::write()
{
    // const word timeDesc =
    //     useTimeName ? time_.timeName() : Foam::name(time_.timeIndex());

    const word timeDesc = "_" +
    (
        printf_.empty()
      ? Foam::name(time_.timeIndex())
      : word::printf(printf_, time_.timeIndex())
    );

    const scalar timeValue = time_.value();

    update();

    if (meshes_.empty() || (!doInternal_ && !doBoundary_))
    {
        // Skip
        return true;
    }


    fileName vtkName = time_.globalCaseName();

    vtk::vtmWriter vtmMultiRegion;

    Info<< name() << " output Time: " << time_.timeName() << nl;

    label regioni = 0;
    for (const word& regionName : meshes_.sortedToc())
    {
        fileName regionPrefix;
        if (regionName != polyMesh::defaultRegion)
        {
            regionPrefix = regionName;
        }

        auto& meshProxy = meshSubsets_[regioni];
        auto& vtuMeshCells = vtuMappings_[regioni];
        ++regioni;

        const fvMesh& baseMesh = meshProxy.baseMesh();

        wordHashSet acceptField(baseMesh.names<void>(selectFields_));

        // Prune restart fields
        acceptField.filterKeys
        (
            [](const word& k){ return k.endsWith("_0"); },
            true // prune
        );

        const label nVolFields =
        (
            (doInternal_ || doBoundary_)
          ? baseMesh.count
            (
                stringListOps::foundOp<word>(fieldTypes::volume),
                acceptField
            )
          : 0
        );

        // Undecided if we want to automatically support DimensionedFields
        // or only on demand:
        const label nDimFields = 0;
        // (
        //     (doInternal_ || doBoundary_)
        //   ? baseMesh.count
        //     (
        //         stringListOps::foundOp<word>(fieldTypes::internal),
        //         acceptField
        //     )
        //   : 0
        // );


        // Setup for the vtm writer.
        // For legacy format, the information added is simply ignored.

        fileName vtmOutputBase
        (
            outputDir_/regionPrefix/vtkName + timeDesc
        );

        // Combined internal + boundary in a vtm file
        vtk::vtmWriter vtmWriter;

        // Collect individual boundaries into a vtm file
        vtk::vtmWriter vtmBoundaries;

        // Setup the internal writer
        autoPtr<vtk::internalWriter> internalWriter;

        // Interpolator for volume and dimensioned fields
        autoPtr<volPointInterpolation> pInterp;

        if (doInternal_)
        {
            if (interpolate_)
            {
                pInterp.reset(new volPointInterpolation(meshProxy.mesh()));
            }

            if (vtuMeshCells.empty())
            {
                // Use the appropriate mesh (baseMesh or subMesh)
                vtuMeshCells.reset(meshProxy.mesh());
            }

            internalWriter = autoPtr<vtk::internalWriter>::New
            (
                meshProxy.mesh(),
                vtuMeshCells,
                writeOpts_,
                // Output name for internal
                (
                    writeOpts_.legacy()
                  ? vtmOutputBase
                  : (vtmOutputBase / "internal")
                ),
                Pstream::parRun()
            );

            Info<< "    Internal  : "
                << time_.relativePath(internalWriter->output())
                << endl;

            // No sub-block for internal
            vtmWriter.append_vtu
            (
                "internal",
                vtmOutputBase.name()/"internal"
            );

            internalWriter->writeTimeValue(timeValue);
            internalWriter->writeGeometry();
        }


        // Setup the patch writers

        const polyBoundaryMesh& patches = meshProxy.mesh().boundaryMesh();

        PtrList<vtk::patchWriter> patchWriters;
        PtrList<PrimitivePatchInterpolation<primitivePatch>> patchInterps;

        labelList patchIds;
        if (doBoundary_)
        {
            patchIds = getSelectedPatches(patches);
        }

        if (oneBoundary_ && patchIds.size())
        {
            auto writer = autoPtr<vtk::patchWriter>::New
            (
                meshProxy.mesh(),
                patchIds,
                writeOpts_,
                // Output name for one patch: "boundary"
                (
                    writeOpts_.legacy()
                  ? (outputDir_/regionPrefix/"boundary"/"boundary" + timeDesc)
                  : (vtmOutputBase / "boundary")
                ),
                Pstream::parRun()
            );

            // No sub-block for one-patch
            vtmWriter.append_vtp
            (
                "boundary",
                vtmOutputBase.name()/"boundary"
            );

            Info<< "    Boundaries: "
                << time_.relativePath(writer->output()) << nl;


            writer->writeTimeValue(timeValue);
            writer->writeGeometry();

            // Transfer writer to list for later use
            patchWriters.resize(1);
            patchWriters.set(0, writer);

            // Avoid patchInterpolation for each sub-patch
            patchInterps.resize(1); // == nullptr
        }
        else if (patchIds.size())
        {
            patchWriters.resize(patchIds.size());
            if (interpolate_)
            {
                patchInterps.resize(patchIds.size());
            }

            label nPatchWriters = 0;
            label nPatchInterps = 0;

            for (const label patchId : patchIds)
            {
                const polyPatch& pp = patches[patchId];

                auto writer = autoPtr<vtk::patchWriter>::New
                (
                    meshProxy.mesh(),
                    labelList(one(), pp.index()),
                    writeOpts_,
                    // Output name for patch: "boundary"/name
                    (
                        writeOpts_.legacy()
                      ?
                        (
                            outputDir_/regionPrefix/pp.name()
                          / (pp.name()) + timeDesc
                        )
                      : (vtmOutputBase / "boundary" / pp.name())
                    ),
                    Pstream::parRun()
                );

                if (!nPatchWriters)
                {
                    vtmWriter.beginBlock("boundary");
                    vtmBoundaries.beginBlock("boundary");
                }

                vtmWriter.append_vtp
                (
                    pp.name(),
                    vtmOutputBase.name()/"boundary"/pp.name()
                );

                vtmBoundaries.append_vtp
                (
                    pp.name(),
                    "boundary"/pp.name()
                );

                Info<< "    Boundary  : "
                    << time_.relativePath(writer->output()) << nl;

                writer->writeTimeValue(timeValue);
                writer->writeGeometry();

                // Transfer writer to list for later use
                patchWriters.set(nPatchWriters++, writer);

                if (patchInterps.size())
                {
                    patchInterps.set
                    (
                        nPatchInterps++,
                        new PrimitivePatchInterpolation<primitivePatch>(pp)
                    );
                }
            }

            if (nPatchWriters)
            {
                vtmWriter.endBlock("boundary");
                vtmBoundaries.endBlock("boundary");
            }

            patchWriters.resize(nPatchWriters);
            patchInterps.resize(nPatchInterps);
        }

        // CellData
        {
            if (internalWriter.valid())
            {
                // cellIds + procIds (parallel)
                internalWriter->beginCellData
                (
                    (writeIds_ ? 1 + (internalWriter->parallel() ? 1 : 0) : 0)
                  + (internalWriter->parallel() ? 1 : 0)
                  + nVolFields + nDimFields
                );

                // Write cellID field + procID (parallel only)
                if (writeIds_)
                {
                    internalWriter->writeCellIDs();
                    internalWriter->writeProcIDs(); // parallel only
                }
            }

            if (nVolFields)
            {
                for (vtk::patchWriter& writer : patchWriters)
                {
                    writer.beginCellData
                    (
                        (writeIds_ ? 1 : 0)
                      + nVolFields
                    );
                    if (writeIds_)
                    {
                        writer.writePatchIDs();
                    }
                }
            }

            writeAllVolFields
            (
                internalWriter,
                patchWriters,
                meshProxy,
                acceptField
            );

            // writeAllDimFields
            // (
            //     internalWriter,
            //     meshProxy,
            //     acceptField
            // );

            // End CellData is implicit
        }


        // PointData
        // - only construct pointMesh on request since it constructs
        //   edge addressing
        if (interpolate_)
        {
            // Begin PointData
            if (internalWriter.valid())
            {
                internalWriter->beginPointData
                (
                    nVolFields + nDimFields
                );
            }

            forAll(patchWriters, writeri)
            {
                const label nPatchFields =
                (
                    writeri < patchInterps.size() && patchInterps.set(writeri)
                  ? nVolFields
                  : 0
                );

                if (nPatchFields)
                {
                    patchWriters[writeri].beginPointData(nPatchFields);
                }
            }

            writeAllVolFields
            (
                internalWriter, pInterp,
                patchWriters,   patchInterps,
                meshProxy,
                acceptField
            );

            // writeAllDimFields
            // (
            //     internalWriter, pInterp,
            //     meshProxy,
            //     acceptField
            // );

            // writeAllPointFields
            // (
            //     internalWriter,
            //     patchWriters,
            //     meshProxy,
            //     acceptField
            // );

            // End PointData is implicit
        }


        // Finish writers
        if (internalWriter.valid())
        {
            internalWriter->close();
        }

        for (vtk::patchWriter& writer : patchWriters)
        {
            writer.close();
        }

        pInterp.clear();
        patchWriters.clear();
        patchInterps.clear();


        // Collective output

        if (Pstream::master())
        {
            // Naming for vtm, file series etc.
            fileName outputName(vtmOutputBase);

            if (writeOpts_.legacy())
            {
                if (doInternal_)
                {
                    // Add to file-series and emit as JSON

                    outputName.ext(vtk::legacy::fileExtension);

                    fileName seriesName(vtk::seriesWriter::base(outputName));

                    vtk::seriesWriter& series = series_(seriesName);

                    // First time?
                    // Load from file, verify against filesystem,
                    // prune time >= currentTime
                    if (series.empty())
                    {
                        series.load(seriesName, true, timeValue);
                    }

                    series.append(timeValue, timeDesc);
                    series.write(seriesName);
                }
            }
            else
            {
                if (vtmWriter.size())
                {
                    // Emit ".vtm"

                    outputName.ext(vtmWriter.ext());

                    vtmWriter.setTime(timeValue);
                    vtmWriter.write(outputName);

                    fileName seriesName(vtk::seriesWriter::base(outputName));

                    vtk::seriesWriter& series = series_(seriesName);

                    // First time?
                    // Load from file, verify against filesystem,
                    // prune time >= currentTime
                    if (series.empty())
                    {
                        series.load(seriesName, true, timeValue);
                    }

                    series.append(timeValue, outputName);
                    series.write(seriesName);

                    // Add to multi-region vtm
                    vtmMultiRegion.add(regionName, regionPrefix, vtmWriter);
                }

                if (vtmBoundaries.size())
                {
                    // Emit "boundary.vtm" with collection of boundaries

                    // Naming for vtm
                    fileName outputName(vtmOutputBase / "boundary");
                    outputName.ext(vtmBoundaries.ext());

                    vtmBoundaries.setTime(timeValue);
                    vtmBoundaries.write(outputName);
                }
            }
        }
    }


    // Emit multi-region vtm
    if (Pstream::master() && meshes_.size() > 1)
    {
        fileName outputName
        (
            outputDir_/vtkName + "-regions" + timeDesc + ".vtm"
        );

        vtmMultiRegion.setTime(timeValue);
        vtmMultiRegion.write(outputName);

        fileName seriesName(vtk::seriesWriter::base(outputName));

        vtk::seriesWriter& series = series_(seriesName);

        // First time?
        // Load from file, verify against filesystem,
        // prune time >= currentTime
        if (series.empty())
        {
            series.load(seriesName, true, timeValue);
        }

        series.append(timeValue, outputName);
        series.write(seriesName);
    }

    return true;
}


bool Foam::functionObjects::vtkWrite::end()
{
    meshSubsets_.clear();
    vtuMappings_.clear();
    meshes_.clear();

    return true;
}


// ************************************************************************* //
