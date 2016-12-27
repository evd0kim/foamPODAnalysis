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

#include "podWriter.H"
#include "coordSet.H"
#include "OFstream.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class Type>
Foam::autoPtr< Foam::podWriter<Type> > Foam::podWriter<Type>::New
(
    const word& writeType
)
{
    typename wordConstructorTable::iterator cstrIter =
        wordConstructorTablePtr_->find(writeType);

    if (cstrIter == wordConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "podWriter::New(const word&)"
        )   << "Unknown write type "
            << writeType << nl << nl
            << "Valid write types : " << endl
            << wordConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<podWriter<Type> >(cstrIter()());
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::fileName Foam::podWriter<Type>::getBaseName
(
 const word& fieldName
) const
{
    fileName fName("raw");
    
    fName += '_' + fieldName;
    
    return fName;
}

template<class Type>
void Foam::podWriter<Type>::writeTable
(
    const List<Type>& values,
    Ostream& os
) const
{
    forAll(values, listItem)
    {
        write(values[listItem], os);
        if (values[listItem] != values.last())
          {
            os<<nl;
          }
    }
}


template<class Type>
void Foam::podWriter<Type>::writeTable
(
    const List<const List<Type>*>& valuesPtrList,
    Ostream& os
) const
{
     forAll(valuesPtrList, listItem)
     {
       const List<Type>& values = *valuesPtrList[listItem];
       forAll(values, i)
         {   
           write(values[i], os);
           if (values[i] != values.last())
             {
               os<<nl;
             }
         }
     }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::podWriter<Type>::podWriter()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::podWriter<Type>::~podWriter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::podWriter<Type>::write
(
    const List<Field<Type> >& valueSets,
    Ostream& os
) const
{
    List<const Field<Type>*> valueSetPtrs(valueSets.size());
    forAll(valueSetPtrs, i)
    {
        valueSetPtrs[i] = &valueSets[i];
    }
    write(valueSetPtrs, os);
}


template<class Type>
Foam::Ostream& Foam::podWriter<Type>::write
(
    const scalar value,
    Ostream& os
) const
{
  return os << value;
}


template<class Type>
template<class VSType>
Foam::Ostream& Foam::podWriter<Type>::writeVS
(
    const VSType& value,
    Ostream& os
) const
{
    for (direction d=0; d<VSType::nComponents; d++)
    {
        os << value.component(d);
    }
    return os;
}

template<class Type>
Foam::Ostream& Foam::podWriter<Type>::write
(
    const vector& value,
    Ostream& os
) const
{
    return writeVS(value, os);
}


template<class Type>
Foam::Ostream& Foam::podWriter<Type>::write
(
    const sphericalTensor& value,
    Ostream& os
) const
{
    return writeVS(value, os);
}


template<class Type>
Foam::Ostream& Foam::podWriter<Type>::write
(
    const symmTensor& value,
    Ostream& os
) const
{
    return writeVS(value, os);
}


template<class Type>
Foam::Ostream& Foam::podWriter<Type>::write
(
    const tensor& value,
    Ostream& os
) const
{
    return writeVS(value, os);
}


// ************************************************************************* //
