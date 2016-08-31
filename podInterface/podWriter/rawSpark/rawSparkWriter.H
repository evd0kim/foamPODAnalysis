/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::rawSparkWriter

Description

SourceFiles
    rawSparkWriter.C

\*---------------------------------------------------------------------------*/

#ifndef rawSparkWriter_H
#define rawSparkWriter_H

#include "podWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class rawSparkWriter Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class rawSparkWriter
:
    public podWriter<Type>
{

public:

    //- Runtime type information
    TypeName("rawSpark");


    // Constructors

        //- Construct null
        rawSparkWriter();


    //- Destructor
    virtual ~rawSparkWriter();


    // Member Functions

        virtual fileName getFileName
        (
         const scalar&,
         const int,
         const word&
        ) const;

        virtual void write
        (
         const List<const Field<Type>*>&,
         Ostream&
        ) const;

        virtual void write
        (
            const List<const Field<Type> >&,
            Ostream&
        ) const;

        virtual void write
        (
            const bool writeTracks,
            const PtrList<coordSet>&,
            const wordList& valueSetNames,
            const List<List<Field<Type> > >&,
            Ostream&
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "rawSparkWriter.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
