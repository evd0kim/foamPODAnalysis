/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2014 OpenFOAM Foundation
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

#include "functionObjectList.H"
#include "podInterface.H"
#include "volFields.H"
#include "dictionary.H"
#include "podInterface.H"
#include "coordSet.H"

#include "globalIndex.H"
#include "argList.H"

#include "IOstreams.H"
#include "OStringStream.H"
#include "IStringStream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(podInterface, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::podInterface::podInterface
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    name_(name),
    obr_(obr),
    globalPODPath_("."),
    active_(true),
    fieldName_("undefined-fieldName"),
    resultName_("undefined-resultName"),
    timeStart_(0),
    timeEnd_(1),
    chunkNumber_(3),
    writePrecision_(9),
    debug_(true),
    times_(),
    scalarFormatterPtr_(NULL),
    vectorFormatterPtr_(NULL)
{
    // Check if the available mesh is an fvMesh, otherwise deactivate
    if (!isA<fvMesh>(obr_))
    {
        active_ = false;
        WarningIn
        (
            "podInterface::podInterface"
            "("
                "const word&, "
                "const objectRegistry&, "
                "const dictionary&, "
                "const bool"
            ")"
        )   << "No fvMesh available, deactivating." << nl
            << endl;
    }

    read(dict);

    const Time& runTime = obr_.time();
    // Make output directory
    globalPODPath_ = Foam::fileName
      (
       Pstream::parRun()
       ? runTime.path()/".."/"postProcessing"/"POD"
       : runTime.path()/"postProcessing"/"POD"
       );

    mkDir(globalPODPath_);
    mkDir(globalPODPath_/"mapping");
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::podInterface::~podInterface()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::podInterface::read(const dictionary& dict)
{
    if (active_)
    {
        dict.lookup("fieldName") >> fieldName_;
        dict.lookup("resultName") >> resultName_;
        dict.lookup("timeStart") >> timeStart_;
        dict.lookup("timeEnd") >> timeEnd_;
        dict.lookup("chunkNumber") >> chunkNumber_;
        dict.lookup("writePrecision") >> writePrecision_;

        debug_ = dict.lookupOrDefault("debug", false);
        
        if (resultName_ == "none")
        {
            resultName_ = "mode" + fieldName_ ;
        }
        // Define the surface formatter
        // Optionally defined extra controls for the output formats
        scalarFormatterPtr_ = podWriter<scalar>::New
          (
           "rawSpark"
           );
        vectorFormatterPtr_ = podWriter<vector>::New
          (
           "rawSpark"
           );
    }

}


void Foam::podInterface::execute()
{

}


void Foam::podInterface::end()
{
     if (active_)
     {
      // saving snapshot time moment
      if (Pstream::master() || !Pstream::parRun())
        {         
          Foam::OFstream mapTimesFile(
                                      globalPODPath_/
                                      "mapping"/
                                      "timeMap.spark"
                                      );
          // no way except for loop
          forAll(times_, idx)
            {
              mapTimesFile<<times_[idx]<<nl;
            }
        }
     }
}


void Foam::podInterface::timeSet()
{
    // Do nothing
}


void Foam::podInterface::write()
{
    if (active_)
    {
      // saving snapshot time moment
      if (Pstream::master() || !Pstream::parRun())
        {
          times_.append(obr_.time().value());
        }
      //-

      // approaching to main task of the library
      if(Pstream::parRun())
        {
          Info<<"Proper Orthogonal Decomposition library. Writing in parallel mode>>"<<nl
              <<"    Field:"<<fieldName_<<nl;

          const scalarField& vf = obr_.lookupObject<scalarField>(fieldName_);        
          int fieldSize = vf.size();
        
          //mapping cells file

          writeMappingFile(Pstream::myProcNo(), fieldSize);
          
          int procChunk = chunkNumber_*(Pstream::myProcNo());

          writeField(procChunk, vf);
          
          if(debug_)
            {
              Pout<<"Chunk field size = " << fieldSize<<nl
                  <<"Proc chunk = "<<procChunk<<nl
                  <<"Proc ID = "<<Pstream::myProcNo()<<endl;
            }
        }
      else
        {
          Info<<"Proper Orthogonal Decomposition library>>"<<nl
              <<"    Field:"<<fieldName_<<nl;

          const scalarField& vf = obr_.lookupObject<scalarField>(fieldName_);        

          int chunkIndex = 0;

          writeField(chunkIndex, vf);
        };
    }                         

}

void Foam::podInterface::writeChunkedField(int& index, const scalarField& chunkedField)
{
  scalar currentTime = obr_.time().value();

  List<const scalarField*> scalarValues(1);

  scalarValues[0] = &chunkedField;

  fileName outputFile
    ( 
     globalPODPath_
     / scalarFormatterPtr_().getFileName
     (
      currentTime,
      index, 
      fieldName_
      )
      );
  
  Foam::OFstream outputStream(outputFile);
  outputStream.precision(writePrecision_);        

  scalarFormatterPtr_().write
    (
     scalarValues,
     outputStream
     );

  Pout<<"    Export data file: "<<outputFile<<nl
      <<"    Output data precision:"<<outputStream.precision()<<nl<<nl;

}

void Foam::podInterface::writeField(int& chunkIndex, const scalarField& field)
{
  int chunkSize = field.size()/chunkNumber_;
  int resultChunk = 0;
  for(int Idx = 0; Idx < chunkNumber_ - 1; Idx++)
    {
      scalarField partToWrite =  
        scalarField::subField(
                              field, 
                              chunkSize, 
                              Idx*chunkSize
                              );
      if(debug_)
        {
          Pout<<"    Local Idx = "<< Idx<<nl
              <<"    Local Size = "<< partToWrite.size()<<nl
              <<"    Global Chunk = "<<resultChunk<<nl;                
        }
      
      resultChunk = Idx + chunkIndex;

      writeChunkedField( resultChunk, partToWrite);
    }

  scalarField partToWrite =  
    scalarField::subField(
                          field, 
                          field.size() % chunkNumber_ + chunkSize, 
                          chunkSize*(chunkNumber_-1)
                          );

  resultChunk = (chunkNumber_ - 1) + chunkIndex;
  writeChunkedField(resultChunk, partToWrite);

  if(debug_)
    {
      Pout<<"    Local Idx = "<<"Last"<<nl                                
          <<"    Local Size = "<< partToWrite.size()<<nl
          <<"    Global Chunk = "<<resultChunk<<nl;                
    }
}

void Foam::podInterface::writeMappingFile
(
 const int& procNumber, 
 const int& fieldSize
 )
{ 
  const fvMesh& mesh = dynamic_cast<const fvMesh&>(obr_);

  std::ostringstream ostr;
  ostr <<"map"<<"Proc"<<procNumber<<".spark";
  std::string converted = ostr.str();

  fileName fName
    (
     globalPODPath_/
     "mapping"/
     converted
     );

  Foam::OFstream mappingFileStream(fName);
  
  //read local cell addressing
  labelIOList localCellProcAddr
    (
     IOobject
     (
      "cellProcAddressing",
      mesh.facesInstance(),
      mesh.meshSubDir,
      mesh,
      IOobject::MUST_READ,
      IOobject::NO_WRITE
      )
     );

  //Pout<<" output"<<localCellProcAddr<<nl;

  for(int i = 0; i < fieldSize; i++)
    {
      //Pout<<"    processor local cellI = "<<i<<nl<<
      //  "    processor global globalCellI = "<<localCellProcAddr[i]<<nl;
      mappingFileStream<< localCellProcAddr[i];
      if (i != fieldSize - 1)
        {
          mappingFileStream<<nl;
        }
    }
}

// ************************************************************************* //
