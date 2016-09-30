/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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

Class
    Foam::podInterface

Group
    grpFVFunctionObjects

Description
    This function object calculates the magnitude of a field.  The operation
    can be applied to any volume or surface fieldsm and the output is a
    volume or surface scalar field.

SourceFiles
    podInterface.C
    IOpodInterface.H

\*---------------------------------------------------------------------------*/

#ifndef podInterface_H
#define podInterface_H

#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "pointFieldFwd.H"
#include "OFstream.H"
#include "Switch.H"
#include "podWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class polyMesh;
class mapPolyMesh;
class dimensionSet;

/*---------------------------------------------------------------------------*\
                           Class podInterface Declaration
\*---------------------------------------------------------------------------*/

class podInterface
{
    // Private data

        //- Name of this podInterface object
        word name_;

        //- Reference to the database
        const objectRegistry& obr_;

        fileName globalPODPath_;

        //- On/off switch
        bool active_;

        //- Name of field to process
        word fieldName_;

        //- Name of result field
        word resultName_;

        //-
        scalar timeStart_;
        
        //-
        scalar timeEnd_;
  //-
        int chunkSize_;

        int chunkNumber_;
  //-
        int writePrecision_;

        bool debug_;
  
        DynamicList<scalar> times_;

//- File podWriter for scalar data
        autoPtr<podWriter<scalar> > scalarFormatterPtr_;

//- File podWriter for vector data
        autoPtr<podWriter<vector> > vectorFormatterPtr_;
                                                  

    // Private Member Functions

        //- Helper function to create/store/return the mag field
        template<class FieldType>
        FieldType& magField(const word& magName, const dimensionSet& dims);

        //- Helper function to calculate the magnitude of different field types
        template<class Type>
        void calc
        (
            const word& fieldName,
            const word& resultName,
            bool& processed
        );

        //- Write sampled fieldName on surface and on outputDir path
        // template<class Type>
        // void writeSurface
        // (
        //     const Field<Type>& values,
        //     const label surfI,
        //     const word& fieldName,
        //     const fileName& outputDir
        // );

        void writeChunkedField
        (
         int& index, 
         const scalarField& chunkedField
         );

        void writeField
        (
         int& chunkIndex, 
         const scalarField& field
         );

        void writeMappingFile
        (
         const int& procNumber, 
         const int& fieldSize
         );

        //- Disallow default bitwise copy construct
        podInterface(const podInterface&);

        //- Disallow default bitwise assignment
        void operator=(const podInterface&);


public:

    //- Runtime type information
    TypeName("podInterface");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        podInterface
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    //- Destructor
    virtual ~podInterface();


    // Member Functions

        //- Return name of the set of podInterface
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the podInterface data
        virtual void read(const dictionary&);

        //- Execute, currently does nothing
        virtual void execute();

        //- Execute at the final time-loop, currently does nothing
        virtual void end();

        //- Called when time was set at the end of the Time::operator++
        virtual void timeSet();

        //- Calculate the podInterface and write
        virtual void write();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&)
        {}

        //- Update for changes of mesh
        virtual void movePoints(const polyMesh&)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "podInterfaceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
