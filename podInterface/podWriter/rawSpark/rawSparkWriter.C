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

\*---------------------------------------------------------------------------*/

#include "rawSparkWriter.H"
#include "coordSet.H"
#include "fileName.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::rawSparkWriter<Type>::rawSparkWriter()
:
    podWriter<Type>()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::rawSparkWriter<Type>::~rawSparkWriter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::fileName Foam::rawSparkWriter<Type>::getFileName
(
 const scalar& time,
 const int chunkNo,
 const word& fieldName
) const
{
  std::ostringstream ostr;
  ostr <<"_"<<time<<"_"<< chunkNo<<".spark";
  std::string converted = ostr.str();

  fileName fName(this->getBaseName(fieldName));
  
  fName += converted;
    
  return fName;
}

template<class Type>
void Foam::rawSparkWriter<Type>::write
(
    const List<const Field<Type>*>& valueSets,
    Ostream& os
) const
{
    // Collect sets into columns
    List<const List<Type> *> columns(valueSets.size());
 
    forAll(valueSets, i)
    {
      columns[i] = valueSets[i];
    }

    this->writeTable(columns, os);
   
}

template<class Type>
void Foam::rawSparkWriter<Type>::write
(
 const List<const Field<Type> >& valueSets,
 Ostream& os
) const
{  

}

template<class Type>
void Foam::rawSparkWriter<Type>::write
(
 const bool writeTracks,
 const PtrList<coordSet>& points,
 const wordList& valueSetNames,
 const List<List<Field<Type> > >& valueSets,
 Ostream& os
) const
{

}

// ************************************************************************* //
