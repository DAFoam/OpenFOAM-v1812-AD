/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

Class
    Foam::functionObjectList

Description
    List of function objects with start(), execute() and end() functions
    that is called for each object.

See also
    Foam::functionObject
    Foam::functionObjects::timeControl

SourceFiles
    functionObjectList.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjectList_H
#define functionObjectList_H

#include "PtrList.H"
#include "functionObject.H"
#include "SHA1Digest.H"
#include "HashTable.H"
#include "IOdictionary.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class argList;
class mapPolyMesh;
class wordRe;

/*---------------------------------------------------------------------------*\
                     Class functionObjectList Declaration
\*---------------------------------------------------------------------------*/

class functionObjectList
:
    private PtrList<functionObject>
{
    // Private data

        //- A list of SHA1 digests for the function object dictionaries
        List<SHA1Digest> digests_;

        //- Quick lookup of the index into functions/digests
        HashTable<label> indices_;

        const Time& time_;

        //- The parent dictionary containing a "functions" entry
        //  This entry can either be a list or a dictionary of
        //  functionObject specifications
        const dictionary& parentDict_;

        //- Function object properties - stores state information
        mutable autoPtr<IOdictionary> stateDictPtr_;

        //- Switch for the execution of the functionObjects
        bool execution_;

        //- Tracks if read() was called while execution is on
        bool updated_;


    // Private Member Functions

        //- Create state dictionary
        void createStateDict() const;

        //- Remove and return the function object pointer by name,
        //- and returns the old index via the parameter.
        //  Returns a nullptr (and index -1) if it didn't exist
        functionObject* remove(const word& key, label& oldIndex);

        //- Search the specified directory for functionObject
        //- configuration files, add to the given map and recurse
        static void listDir(const fileName& dir, wordHashSet& available);

        //- No copy construct
        functionObjectList(const functionObjectList&) = delete;

        //- No copy assignment
        void operator=(const functionObjectList&) = delete;


public:

    // Static data members

        //- Default relative path ("caseDicts/postProcessing") to the
        //- directory structure containing functionObject dictionary files.
        static fileName functionObjectDictPath;


    // Constructors

        //- Construct from Time and the execution setting.
        //  The functionObject specifications are read from the controlDict
        functionObjectList
        (
            const Time& runTime,
            const bool execution=true
        );

        //- Construct from Time, a dictionary with "functions" entry
        //  and the execution setting.
        //  \param[in]  runTime - the other Time instance to construct from
        //  \param[in]  parentDict - the parent dictionary containing
        //    a "functions" entry, which can either be a list or a dictionary
        //    of functionObject specifications.
        //  \param[in]  execution - whether the function objects should execute
        //    or not. Default: true.
        functionObjectList
        (
            const Time& runTime,
            const dictionary& parentDict,
            const bool execution=true
        );

        //- Construct and return a functionObjectList for an application.
        //  If the "dict" argument is specified the functionObjectList is
        //  constructed from that dictionary which is returned as
        //  controlDict otherwise the functionObjectList is constructed
        //  from the "functions" sub-dictionary of "system/controlDict"
        static autoPtr<functionObjectList> New
        (
            const argList& args,
            const Time& runTime,
            dictionary& controlDict,
            HashSet<wordRe>& requiredFields
        );


    //- Destructor
    ~functionObjectList() = default;


    // Member Functions

        //- Return the number of elements in the List.
        using PtrList<functionObject>::size;

        //- Return true if the List is empty (ie, size() is zero).
        using PtrList<functionObject>::empty;

        //- Access to the functionObjects
        using PtrList<functionObject>::operator[];

        //- Reset/read state dictionary for current time
        void resetState();

        //- Return the state dictionary
        IOdictionary& stateDict();

        //- Return const access to the state dictionary
        const IOdictionary& stateDict() const;

        //- Clear the list of function objects
        void clear();

        //- Find the ID of a given function object by name
        label findObjectID(const word& name) const;

        //- Print a list of functionObject configuration files in the
        //- directories located using
        //- Foam::findEtcDirs("caseDicts/postProcessing")
        //
        //  -# \b user settings
        //    - ~/.OpenFOAM/{PROJECT_API}/"caseDicts/postProcessing"
        //    - ~/.OpenFOAM/"caseDicts/postProcessing"
        //  -# \b group settings
        //    - $WM_PROJECT_SITE/{PROJECT_API}/"etc/caseDicts/postProcessing"
        //    - $WM_PROJECT_SITE/"etc/caseDicts/postProcessing"
        //  -# \b other (shipped) settings
        //    - $WM_PROJECT_DIR/etc/"caseDicts/postProcessing"
        //
        // Where {PROJECT_API} is the value of the OPENFOAM define.
        // See further notes in Foam::findEtcEntries()
        static void list();

        //- Find a functionObject dictionary file in the case
        //- \<system\> directory or any directory located using
        //- Foam::findEtcDirs("caseDicts/postProcessing")
        //
        //  \return The path of the functionObject dictionary file found
        //     or an empty path
        static fileName findDict(const word& funcName);

        //- Read the specified functionObject configuration dictionary parsing
        //- the optional arguments included in the name 'funcNameArgs0',
        //- inserting 'field' or 'fields' entries as required and merging the
        //- resulting functionObject dictionary into 'functionsDict'.  Any
        //- fields required to execute the functionObject are added to
        //- 'requiredFields'
        //
        //  Uses functionObjectList::findDict() for searching
        static bool readFunctionObject
        (
            const string& funcNameArgs0,
            dictionary& functionsDict,
            HashSet<wordRe>& requiredFields,
            const word& region = word::null
        );

        //- Read and set the function objects if their data have changed
        bool read();

        //- Switch the function objects on
        void on();

        //- Switch the function objects off
        void off();

        //- Return the execution status (on/off) of the function objects
        bool status() const;

        //- Called at the start of the time-loop
        bool start();

        //- Called at each ++ or += of the time-loop.
        //  postProcess overrides the usual executeControl behaviour and
        //  forces execution (used in post-processing mode)
        bool execute();

        //- Execute function objects using the specified subIndex.
        //  \param subIndex an execution sub-index corresponding to a
        //      sub-cycle or something similar
        bool execute(const label subIndex);

        //- Execute a subset of function objects using the specified subIndex.
        //  \param functionNames names or regex of existing functions to
        //      execute
        //  \param subIndex an execution sub-index corresponding to a
        //      sub-cycle or something similar
        bool execute(const UList<wordRe>& functionNames, const label subIndex);

        //- Called when Time::run() determines that the time-loop exits
        bool end();

        //- Called at the end of Time::adjustDeltaT() if adjustTime is true
        bool adjustTimeStep();

        //- Did any file get changed during execution?
        bool filesModified() const;

        //- Update for changes of mesh
        void updateMesh(const mapPolyMesh& mpm);

        //- Update for changes of mesh
        void movePoints(const polyMesh& mesh);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
