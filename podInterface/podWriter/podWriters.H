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

InClass
    Foam::podWriter

SourceFiles
    podWriters.C

\*---------------------------------------------------------------------------*/

#ifndef podWriters_H
#define podWriters_H

#include "podWriter.H"
#include "fieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Only used internally
#define makeTypeSetWritersTypeName(typePodWriter, dataType)                      \
                                                                              \
    defineNamedTemplateTypeNameAndDebug(typePodWriter< dataType >, 0)

// Sometimes used externally
#define makeSetWritersTypeName(typePodWriter)                                    \
                                                                              \
    makeTypeSetWritersTypeName(typePodWriter, scalar);                           \
    makeTypeSetWritersTypeName(typePodWriter, vector);                           \
    makeTypeSetWritersTypeName(typePodWriter, sphericalTensor);                  \
    makeTypeSetWritersTypeName(typePodWriter, symmTensor);                       \
    makeTypeSetWritersTypeName(typePodWriter, tensor)

// Define type info for single dataType template instantiation (eg, vector)
#define makeSetWriterType(typePodWriter, dataType)                               \
                                                                              \
    defineNamedTemplateTypeNameAndDebug(typePodWriter< dataType >, 0);           \
    addTemplatedToRunTimeSelectionTable                                       \
    (                                                                         \
        podWriter, typePodWriter, dataType, word                                    \
    )


// Define type info for scalar, vector etc. instantiations
#define makeSetWriters(typePodWriter)                                            \
                                                                              \
    makeSetWriterType(typePodWriter, scalar);                                    \
    makeSetWriterType(typePodWriter, vector);                                    \
    makeSetWriterType(typePodWriter, sphericalTensor);                           \
    makeSetWriterType(typePodWriter, symmTensor);                                \
    makeSetWriterType(typePodWriter, tensor)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
