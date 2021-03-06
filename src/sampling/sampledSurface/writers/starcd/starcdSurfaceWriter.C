/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2018 OpenCFD Ltd.
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

#include "starcdSurfaceWriter.H"
#include "MeshedSurfaceProxy.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "makeSurfaceWriterMethods.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    makeSurfaceWriterType(starcdSurfaceWriter);
}

// Field writing implementation
#include "starcdSurfaceWriterImpl.C"

// Field writing methods
defineSurfaceWriterWriteFields(Foam::starcdSurfaceWriter);


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::fileName Foam::starcdSurfaceWriter::write
(
    const fileName& outputDir,
    const fileName& surfaceName,
    const meshedSurf& surf,
    const bool verbose
) const
{
    // geometry:  rootdir/time/surfaceName.{raw,vrt,inp}

    fileName outputFile(outputDir/surfaceName + ".inp");

    if (verbose)
    {
        Info<< "Writing geometry to " << outputFile << endl;
    }

    if (!isDir(outputFile.path()))
    {
        mkDir(outputFile.path());
    }

    MeshedSurfaceProxy<face>(surf.points(), surf.faces()).write
    (
        outputFile
    );

    return outputFile;
}


// ************************************************************************* //
