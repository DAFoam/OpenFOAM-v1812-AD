/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "fileStat.H"
#include "IOstreams.H"
#include "timer.H"

#include <signal.h>
#include <unistd.h>
#ifndef darwin
    #include <sys/sysmacros.h>
#endif

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileStat::fileStat()
:
    isValid_(false)
{}


Foam::fileStat::fileStat
(
    const char* fName,
    const bool followLink,
    const unsigned int maxTime
)
:
    isValid_(false)
{
    if (!fName || !fName[0])
    {
        return;
    }

    // Work on volatile
    volatile bool locIsValid = false;

    timer myTimer(maxTime);

    if (!timedOut(myTimer))
    {
        if (followLink)
        {
            locIsValid = (::stat(fName, &status_) == 0);
        }
        else
        {
            locIsValid = (::lstat(fName, &status_) == 0);
        }
    }

    // Copy into (non-volatile, possible register based) member var
    isValid_ = locIsValid;
}


Foam::fileStat::fileStat
(
    const fileName& fName,
    const bool followLink,
    const unsigned int maxTime
)
:
    fileStat(fName.c_str(), followLink, maxTime)
{}


Foam::fileStat::fileStat(Istream& is)
{
    is >> *this;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::fileStat::size() const
{
    return isValid_ ? label(status_.st_size) : 0;
}


time_t Foam::fileStat::modTime() const
{
    return isValid_ ? status_.st_mtime : 0;
}


double Foam::fileStat::dmodTime() const
{
    return
    (
        isValid_
      ?
        #ifdef darwin
        (status_.st_mtime + 1e-9*status_.st_mtimespec.tv_nsec)
        #else
        (status_.st_mtime + 1e-9*status_.st_mtim.tv_nsec)
        #endif
      : 0
    );
}


bool Foam::fileStat::sameDevice(const fileStat& other) const
{
    return
        isValid_
     && (
            major(status_.st_dev) == major(other.status_.st_dev)
         && minor(status_.st_dev) == minor(other.status_.st_dev)
        );
}


bool Foam::fileStat::sameINode(const fileStat& other) const
{
    return isValid_ && (status_.st_ino == other.status_.st_ino);
}


bool Foam::fileStat::sameINode(const label iNode) const
{
    return isValid_ && (status_.st_ino == ino_t(iNode));
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, fileStat& fStat)
{
    FixedList<label, 13> stat(is);

    fStat.isValid_ = stat[0];

    dev_t st_dev = makedev(stat[1], stat[2]);
    fStat.status_.st_dev = st_dev;

    fStat.status_.st_ino = stat[3];
    fStat.status_.st_mode = stat[4];
    fStat.status_.st_uid = stat[5];
    fStat.status_.st_gid = stat[6];

    dev_t st_rdev = makedev(stat[7], stat[8]);
    fStat.status_.st_rdev = st_rdev;

    fStat.status_.st_size = stat[9];
    fStat.status_.st_atime = stat[10];
    fStat.status_.st_mtime = stat[11];
    fStat.status_.st_ctime = stat[12];

    is.check(FUNCTION_NAME);
    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const fileStat& fStat)
{
    FixedList<label, 13> stat;

    stat[0] = label(fStat.isValid_);
    stat[1] = label(major(fStat.status_.st_dev));
    stat[2] = label(minor(fStat.status_.st_dev));
    stat[3] = label(fStat.status_.st_ino);
    stat[4] = label(fStat.status_.st_mode);
    stat[5] = label(fStat.status_.st_uid);
    stat[6] = label(fStat.status_.st_gid);
    stat[7] = label(major(fStat.status_.st_rdev));
    stat[8] = label(minor(fStat.status_.st_rdev));
    stat[9] = label(fStat.status_.st_size);
    stat[10] = label(fStat.status_.st_atime);
    stat[11] = label(fStat.status_.st_mtime);
    stat[12] = label(fStat.status_.st_ctime);

    return os << stat;
}


// ************************************************************************* //
