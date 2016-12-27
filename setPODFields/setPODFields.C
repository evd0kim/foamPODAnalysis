/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

Description
    Set values on a selected set of cells/patchfaces through a dictionary.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "fvCFD.H"
#include "topoSetSource.H"
#include "cellSet.H"
#include "faceSet.H"
#include "volFields.H"

#include "IOstreams.H"
#include "IFstream.H"
#include "OStringStream.H"
#include "IStringStream.H"
#include "OSspecific.H"

#include "SortableList.H"
#include "ListOps.H"

#include "transformGeometricField.H"

using namespace Foam;

/*
template<>
Foam::label readValue(const string& line)
  {
    return readLabel(IStringStream(line)());
  }

template<>
Foam::scalar readValue(const string& line)
{
  return readScalar(IStringStream(line)());
}

template<class Type>
Type readValue(const string& splitted)
{
        Type result;
        
        result =
          readScalar(IStringStream(line)());

        return result;
}

template<class Type>
  DynamicList< Type > readRawFile
(
    const fileName& sourceFile
)
{
  Foam::IFstream inputStream(sourceFile);
              
  Info<<"    Trying to read file: "<<nl
      <<"    "<<inputStream.name()<<endl;

  if (!inputStream.good())
    {
      FatalErrorIn
        (
         "setPODFields::readMode"
         )
        << "Cannot read file " << inputStream.name()
        << exit(FatalError);
    }
  else{
    Info<<"    OK. Read values"<<endl;
    inputStream.defaultPrecision(16);
  }

  DynamicList< Type > values;
              
  while (inputStream.good())
    {
      string line;
      inputStream.getLine(line);
      
      Type x = readValue(IStringStream(line)());
      
      values.append(x);     
    }

  return values;
  }
*/

DynamicList< scalar > readRawFile
(
 const fileName& sourceFile
)
{
  Foam::IFstream inputStream(sourceFile);

  Info<<"    Trying to read file: "<<nl
      <<"    "<<inputStream.name()<<endl;

  if (!inputStream.good())
    {
      FatalErrorIn
        (
         "setPODFields::readMode"
         )
        << "Cannot read file " << inputStream.name()
        << exit(FatalError);
    }
  else{
    Info<<"    OK. Read values"
        <<endl;
    inputStream.defaultPrecision(16);
  }

  DynamicList< scalar > values;
              
  while (inputStream.good())
    {
      string line;
      inputStream.getLine(line);

      wordRe numeric ("[:digit:]*(\.|\,)*[:digit:]*");
      numeric.compile(wordRe::REGEXP);

      if(numeric.match(line))
        values.append(     
                      Foam::readScalar(IStringStream(line)())
                           );
      else
        continue;
      
    }
  
  return values;
}

DynamicList< label > readMappingFile
(
 const fileName& sourceFile
)
{
  Foam::IFstream inputStream(sourceFile);

  Info<<"    Trying to read file: "<<nl
      <<"    "<<inputStream.name()<<endl;

  if (!inputStream.good())
    {
      FatalErrorIn
        (
         "setPODFields::readMode"
         )
        << "Cannot read file " << inputStream.name()
        << exit(FatalError);
    }
  else{
    Info<<"    OK. Read values"
        <<nl<<endl;
    inputStream.defaultPrecision(16);
  }

  DynamicList< label > values;
              
  while (inputStream.good())
    {
      string line;
      inputStream.getLine(line);
      
      if (line=="\n") continue;
      
      label value = Foam::readLabel(IStringStream(line)());
      
      values.append(value);     
    }
  
  return values;
}

DynamicList< scalar > readLineAndSplit
(
 IFstream& is,
 char separator
)
{ 
  DynamicList< scalar > returnList(0);

  if (is.good())
   {
     string line;
     is.getLine(line);
     
     label n = 0;
     std::size_t pos = 0;
     DynamicList<string> splitted;

     std::size_t nPos = 0;
     
     while ( pos != std::string::npos )
       {
         bool found = false;
         while (!found)
           {
             nPos = line.find(separator, pos);
             
             if ((nPos != std::string::npos) && (nPos - pos == 0))
               {
                 pos = nPos + 1;
               }
             else
               {
                 found = true;
               }
           }
         
         nPos = line.find(separator, pos);
         
         if (nPos == std::string::npos)
           {
             splitted.append(line.substr(pos));
             pos = nPos;
             n++;
           }
         else
           {
             splitted.append(line.substr(pos, nPos - pos));
             pos = nPos + 1;
             n++;
           }
       }
  
     forAll(splitted, idx)
       {
         scalar value = Foam::readScalar(IStringStream(splitted[idx])());
         returnList.append(value);     
       }
   }
  
  return returnList;
}

template<class Type>
bool reconstructCellFieldType
(
    const word& fieldTypeDesc,
    const fvMesh& mesh,
    Istream& fieldValueStream,
    const IOdictionary& settingsDict
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    if (fieldTypeDesc != fieldType::typeName)
    {
        return false;
    }

    // reconstruction field name, read from stream
    word fieldName(fieldValueStream);
    
    // reconstruction mode
    word mode = 
      settingsDict.lookupOrDefault<word>("mode", "serial");

    // relative time switch. relative means that goal time is taken between first and last snapshots
    word relativeFlag = 
      settingsDict.lookupOrDefault<word>("relativeTime", "true");

    // the time we are reconstructing for
    scalar goalTime = 
      settingsDict.lookupOrDefault("goalTime", 0.0);

    // amount of modes we take for reconstruction   
    int modesNumber = 
      settingsDict.lookupOrDefault("usingModes", 1);

    // save data  switch. Manage to save modes in each own separate field
    word saveFlag = 
      settingsDict.lookupOrDefault<word>("saveModes", "false");

    // Input for POD -- related files search
    word timeFileName = 
      settingsDict.lookupOrDefault<word>("timeFile", "timeCoeffs.txt");

    // 
    word modesFileName = 
      settingsDict.lookupOrDefault<word>("modesFile", "mode*");
    wordRe modesMaskName (modesFileName);

    if (wordRe::isPattern(modesFileName))
      {
	modesMaskName.compile(wordRe::REGEXP);
      }
    else
      {
	modesMaskName.compile(wordRe::LITERAL);
      }
    
    // Check the current time directory
    
    IOobject fieldHeader
    (
        fieldName,
        mesh.time().timeName(),
        mesh,
        IOobject::MUST_READ
     );

    // Check the "constant" directory
    if (!fieldHeader.headerOk())
    {
      fieldHeader = IOobject
        (
	 fieldName,
	 mesh.time().constant(),
	 mesh,
	 IOobject::MUST_READ
	 );
    }
    
    // Check field exists
    if (fieldHeader.headerOk())
    {
        Info<< "    Reconstructing internal values of "
            << fieldHeader.headerClassName()
            << " " << fieldName << endl;

        fieldType field(fieldHeader, mesh);

        const Type& value = pTraits<Type>(fieldValueStream);
        	
	forAll(field, celli)
        {
	  field[celli] = value;
        }
       
        forAll(field.boundaryField(), patchi)
        {
            field.boundaryField()[patchi] =
                field.boundaryField()[patchi].patchInternalField();
        }

        if (!field.write())
        {
            FatalErrorIn
            (
                "void reconstructCellFieldType"
                "(const fvMesh& mesh, const labelList& selectedCells,"
                "Istream& fieldValueStream)"
            ) << "Failed writing field " << fieldName << endl;
        }
    }
    else
    {
        WarningIn
        (
            "void reconstructCellFieldType"
            "(const fvMesh& mesh, const labelList& selectedCells,"
            "Istream& fieldValueStream)"
        ) << "Field " << fieldName << " not found" << endl;

        // Consume value
        (void)pTraits<Type>(fieldValueStream);
    }

    return true;
}

class reconstructCellField
{

public:

    reconstructCellField()
    {}

    autoPtr<reconstructCellField> clone() const
    {
        return autoPtr<reconstructCellField>(new reconstructCellField());
    }

    class iNew
    {
      const fvMesh& mesh_;
      const IOdictionary& dict_;

    public:

      iNew(const fvMesh& mesh, const IOdictionary& dict)
        :
	mesh_(mesh),
	dict_(dict)
      {}

        autoPtr<reconstructCellField> operator()(Istream& fieldValues) const
        {
            word fieldType(fieldValues);

            if
            (
	     !(
	       reconstructCellFieldType<scalar>
	       (fieldType, mesh_, fieldValues, dict_)
	       || reconstructCellFieldType<vector>
	       (fieldType, mesh_, fieldValues, dict_)
	       )
	     )
	      {
                WarningIn("reconstructCellField::iNew::operator()(Istream& is)")
		  << "field type " << fieldType << " not currently supported"
		  << endl;
	      }

            return autoPtr<reconstructCellField>(new reconstructCellField());
        }
    };
};
  
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
  //#include "createNamedMesh.H"

    Info<< "Reading setPODFieldsDict\n" << endl;

    //Create and initialize dictionary

    IOdictionary setPODFieldsDict
    (
        IOobject
        (
            "setPODFieldsDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    if (setPODFieldsDict.found("reconstructFields"))
      {
        Info<< "Reconstructing field values" << endl;       
	wordReList fieldNames;
	setPODFieldsDict.lookup("reconstructFields") >> fieldNames;
      }

    /*
    //read default values
    if (setPODFieldsDict.found("defaultFieldValues"))
    {
        Info<< "Setting field default values" << endl;
        PtrList<setCellField> defaultFieldValues
        (
            setPODFieldsDict.lookup("defaultFieldValues"),
            setCellField::iNew(mesh, labelList(mesh.nCells()))
        );
        Info<< endl;
    }
    //
    */
    // reconstruction field name  
    word fieldName = setPODFieldsDict.lookupOrDefault<word>("field", "null");
    
    // reconstruction mode
    word mode = setPODFieldsDict.lookupOrDefault<word>("mode", "serial");

    // relative time switch. relative means that goal time is taken between first and last snapshots
    word relativeFlag = setPODFieldsDict.lookupOrDefault<word>("relativeTime", "true");

    // the time we are reconstructing for
    scalar goalTime = setPODFieldsDict.lookupOrDefault("goalTime", 0.0);

    // amount of modes we take for reconstruction   
    int modesNumber = setPODFieldsDict.lookupOrDefault("usingModes", 0);

    // save data  switch. Manage to save modes in each own separate field
    word saveFlag = setPODFieldsDict.lookupOrDefault<word>("saveModes", "false");

    // Input for POD -- related files search
    word timeFileName = setPODFieldsDict.lookupOrDefault<word>("timeFile", "timeCoeffs.txt");

    // 
    word modesFileName = setPODFieldsDict.lookupOrDefault<word>("modesFile", "mode*");
    wordRe modesMaskName (modesFileName);

    if (wordRe::isPattern(modesFileName))
      {
        modesMaskName.compile(wordRe::REGEXP);
      }
    else
      {
        modesMaskName.compile(wordRe::LITERAL);
      }
    
    //int modesNumber(readInt(setPODFieldsDict.lookup("usingModes")));

    Info<< "Setting field values using Proper Orthogonal Decomposition data" << endl;
    
    volScalarField fieldFromTime
      (
        IOobject
        (
            fieldName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
     );

    Foam::scalarField& scalarField = fieldFromTime.internalField();

    fileName inputModesDir
    (
     runTime.path()/"postProcessing"/"POD"/"modes"
     );

    fileName inputMappingDir
      (
       runTime.path()/"postProcessing"/"POD"/"mapping"
     );

    fileName timeMappingFile
      (
       runTime.path()/"postProcessing"/"POD"/"mapping"/"timeMapping"
       );

    fileName timeCoeffsFile
        (
        inputModesDir/timeFileName
        );
   
    // Pre-checks
    if ( modesNumber == 0
         || !isDir(inputModesDir)
         || ( relativeFlag == "false" && !isFile(timeMappingFile) ) 
         || !isFile(timeCoeffsFile)
        )
    {
      if(modesNumber == 0)
        {
          Info()<<nl<<"    The desired reconstruction modes nuber is zero!"
                <<endl;
        }
      else if (!isDir(inputModesDir))
        {
          Info()<<nl<<"    The directory with POD modes is not found!"
                <<endl;
        }
      else if ( relativeFlag == "false" && !isFile(timeMappingFile) )
        {
          Info()<<nl<<"    The absolute time mapping file is not found!"
                <<endl;
        }
        else if ( !isFile(timeCoeffsFile) )
        {
          Info()<<nl<<"    The time coefficients file is not found!"
                <<endl;

        }
      else
        {
          Info()<<nl<<"Error!"<<endl;
        }
    }
    else //everything is OK
        {        
          // Read directory entries into a list. Order makes sense

          Foam::SortableList<fileName> modeEntries
            (
             readDir(inputModesDir, fileName::FILE)
             );

          Info<<"    Found "<<modeEntries.size()<<" files in modes directory"<<endl;

          // Create mapping list
          DynamicList< label > cellsMappingList;

          if(mode == "parallel")
            {
              // Read mapping entries
              // The mapping sequence should save order of processors

              Foam::SortableList<fileName> mapEntries
                ( 
                 readDir(inputMappingDir, fileName::FILE) 
                  );
              
              Info<<"    Proceed in parallel mode. "
                  <<"Found "<<mapEntries.size()<<" files in mapping directory"
                  <<nl<<endl;
                           
              // Gathering cell mapping info from files
              // set up mask for mapping files
              
              wordRe mapFileMask ("mapProc.*");
              mapFileMask.compile(wordRe::REGEXP);
              
              forAll(mapEntries, i)
                {
                  word tmpFileName(mapEntries[i].name(1));
                  // Processor mapping and time mapping are in the same dir so elaborate them
                  // in the same part of code
                  if(mapFileMask.match(tmpFileName))
                    {
                      DynamicList< label > tmpList;
                      tmpList = readMappingFile(
                                                inputMappingDir
                                                /
                                                mapEntries[i]
                                                );
                      cellsMappingList.append(tmpList);
                    }
                  else if (relativeFlag == "false")
                    {
                      if(timeMappingFile == fileName(tmpFileName))
                        {
                          DynamicList< scalar > snapshotTimes;
                          
                          snapshotTimes = readRawFile(timeMappingFile);
                        }
                    }
                }
            }
          else if (mode == "serial")
            {
              for (int i = 0; i<mesh.nCells(); i++)
                {
                  cellsMappingList.append(i);
                }
            }
          else
            {
              FatalErrorIn
                (
                 "setPODFields::int main"
                 )
                <<"Unknown mode for reconstruction. Use serial or parallel mode "<<nl
                << exit(FatalError);
            }

          // Final check for mapping file and interpolated field
          if (cellsMappingList.size()!=scalarField.size())
            {
              FatalErrorIn
                (
                 "setPODFields::int main"
                 )
                <<"The cells mapping list size is not equal to field size "<<nl
                <<"Mapping list size "<<cellsMappingList.size()<<nl
                <<"Field size "<<scalarField.size()<<nl
                << exit(FatalError);
            }
          else{
            Info<<"    Cells mapping: OK"
                <<nl<<endl;
          }
          //-

          // Working with modes
          List<DynamicList< scalar > > modesList(modesNumber);
          List<DynamicList< scalar > > timeCoeffsList(modesNumber);
          int foundModes(0);

          forAll(modeEntries, i)
            {
              word tmpFileName(modeEntries[i].name(1));

              if(modesMaskName.match( modeEntries[i].name() ) && foundModes<modesNumber)
                {                  
                  DynamicList< scalar > values;
                  
                  values = readRawFile(
                                       inputModesDir
                                       /
                                       modeEntries[i]
                                       );
                  
                  if(values.size()==scalarField.size())
                    {                          
                      Info<<"    Mode values: OK"
                          <<"    Adding mode "<<i<<nl<<endl;
                     
                      modesList[foundModes] =  values;
                    }
                  else
                    {
                      FatalErrorIn
                        (
                         "setPODFields::int main"
                         )
                        <<"    Problematic source file >> "<<nl
                        <<"         File size "<<values.size()<<nl
                        <<"         Field size "<<scalarField.size()<<nl
                        << exit(FatalError);
                    }
                  foundModes++;
                }
              else if (modeEntries[i].name() == timeCoeffsFile.name()) 
                {
                  Info<<nl<<"    Time coefficents file has been found. "
                      <<"Take first " << timeCoeffsList.size()<<" modes "<<nl<<endl;

                  char separator = '\t';

                  Foam::IFstream timeCoeffStream( inputModesDir
                                                  /
                                                  modeEntries[i].name()
                                                  );
                  
                  forAll(timeCoeffsList, modeIndex)
                    {
                      timeCoeffsList[modeIndex] = readLineAndSplit(timeCoeffStream, separator);
                      //Info<<"Values "<<timeCoeffsList[modeIndex]<<endl;
                    }
                }                
              else
                {
                }
            }
          
          // saving modes of interest into separated fields to vizualize them all
          if (saveFlag == "true")
            {
              forAll(modesList, modeIndex)
                {
                  std::ostringstream ostr;
                  ostr <<"mode"<<modeIndex;
                  std::string converted = ostr.str();
                  
                  word modeName(converted);
                  
                  Foam::scalarField tmpField( modesList[modeIndex] );

                  volScalarField modeField
                    (
                     IOobject
                     (
                      modeName,
                      runTime.timeName(),
                      mesh,
                      IOobject::NO_READ,
                      IOobject::AUTO_WRITE
                      ),
                     mesh,
                     dimless,                  
                     zeroGradientFvPatchScalarField::typeName
                     );

                  Foam::scalarField& modeScalars = modeField.internalField();                  
                                    
                  forAll(tmpField, idx)
                    {
                      label cellIndex = cellsMappingList[idx];
                      modeScalars[cellIndex] = tmpField[idx];
                    };                

                  modeField.write();
                }
            }
          //-

          //Info<<"    List time coeffs = "
          //    <<timeCoeffsList<<endl;
	  forAll(scalarField, idx)
	  {
	    scalarField[idx] = 0;
	    };

          // interpolating using relative or absolute time
          if (relativeFlag == "true")
            {
              int timePoint(goalTime);
              forAll(modesList, modeIndex)
                {
                  for (int i = 0; i<scalarField.size(); i++)
                    {
                      label cellIndex = cellsMappingList[i];
                      scalarField[cellIndex] += 
                        timeCoeffsList[modeIndex][timePoint]
                        *modesList[modeIndex][i] ;
                    };

                  int cellIndex = cellsMappingList[100];
                  Info()<<"Reconstruction time coeff = "
                        <<timeCoeffsList[modeIndex][timePoint]<<nl
                        <<"Mode value = "
                        <<modesList[modeIndex][100]<<nl
                        <<"Pressure value at point = "
                        <<scalarField[cellIndex]<<nl
                    	<<"Max reconstructed val = "
                    	<<max(scalarField)<<endl;           
                }
            }
          else
            {
            }                 
          // -

          fieldFromTime.write();      
        };

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
