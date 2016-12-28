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
#include "zeroGradientFvPatchFields.H"

using namespace Foam;

typedef struct settingsStruct
{
  word mode;
  word vectorMode;
  word relativeFlag;
  scalar goalTime; 
  int modesNumber;
  word saveFlag; 
  word timeFileName;
  word modesFileName;
};

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

void readPODSettings
(
 settingsStruct& PODSettings, 
 const IOdictionary& settingsDict
 )
{
    // reconstruction mode
    PODSettings.mode = 
      settingsDict.lookupOrDefault<word>("mode", "serial");
    // option for vector reconstruction
    PODSettings.vectorMode = 
      settingsDict.lookupOrDefault<word>("vectorMode", "components");
    // relative time switch. relative means that goal time is taken between first and last snapshots
    PODSettings.relativeFlag = 
      settingsDict.lookupOrDefault<word>("relativeTime", "true");
    // the time we are reconstructing for
    PODSettings.goalTime = 
      settingsDict.lookupOrDefault("goalTime", 0.0);
    // amount of modes we take for reconstruction   
    PODSettings.modesNumber = 
      settingsDict.lookupOrDefault("usingModes", 1);
    // flag to save modes as fields   
    PODSettings.saveFlag = 
      settingsDict.lookupOrDefault<word>("saveModes", "false");
    // Input for POD -- related files search
    PODSettings.timeFileName = 
      settingsDict.lookupOrDefault<word>("timeFile", "timeCoeffs.txt");
    // 
    PODSettings.modesFileName = 
      settingsDict.lookupOrDefault<word>("modesFile", "mode*");
}

void readTimeCoefficients
( 
 const fileName& timeCoeffsFile,
 List<DynamicList< scalar > >& timeCoeffsList,
 char separator
)
{
  Foam::IFstream timeCoeffStream
    (
     timeCoeffsFile
      );

  if ( !isFile(timeCoeffsFile) )
    {
      FatalErrorIn
        (
         "setPODFields::void readTimeCoefficients"
         )
        <<"    The time coefficients file is not found!"<<nl
        <<"    file:"<<timeCoeffsFile<<nl
        << exit(FatalError);
    }

  Info<<"    Time coefficents file has been found. "
      <<"Take first " << timeCoeffsList.size()
      <<" modes."<<endl;
      
  forAll(timeCoeffsList, modeIndex)
    {
      timeCoeffsList[modeIndex] = 
	readLineAndSplit(timeCoeffStream, 
			 separator);
    };

  Info<<"    Timesteps summary "<<timeCoeffsList[0].size()
      <<"."<<endl;
}

template<class Type>
bool reconstructCellFieldType
(
    const word& fieldTypeDesc,
    const fvMesh& mesh,
    const labelList& cellsMappingList,
    Istream& fieldValueStream,
    const settingsStruct& PODSettings,
    const fileName& caseDir
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    if (fieldTypeDesc != fieldType::typeName)
    {
        return false;
    }

    // reconstruction field name, read from stream
    word fieldName(fieldValueStream);

    word mode = PODSettings.mode;
    word relativeFlag = PODSettings.relativeFlag;   
    scalar goalTime = PODSettings.goalTime;
    int modesNumber = PODSettings.modesNumber;
    word saveFlag = PODSettings.saveFlag;
    word timeFileName = PODSettings.timeFileName;
    word modesFileName = PODSettings.modesFileName;

    wordRe modesMaskName (modesFileName);

    if (wordRe::isPattern(modesFileName))
      {
	modesMaskName.compile(wordRe::REGEXP);
      }
    else
      {
	modesMaskName.compile(wordRe::LITERAL);
      }
   
    //Set up file names
    fileName inputModesDir
      (
       caseDir/"postProcessing"/"POD"/"modes"/fieldName
       );

    fileName inputMappingDir
      (
       caseDir/"postProcessing"/"POD"/"mapping"
       );

    fileName timeMappingFile
      (
       caseDir/"postProcessing"/"POD"/"mapping"/"timeMapping"
       );

    fileName timeCoeffsFile
      (
       inputModesDir/timeFileName
       );
    //  
    
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
      Info<<nl
          << "    Reconstructing internal values of "
          << fieldHeader.headerClassName()
          << " " << fieldName << endl;

        fieldType field(fieldHeader, mesh);
	forAll(field, celli)
	  {
	    field[celli] *= 0.0;
	  }

	// Final check for mapping file and interpolated field
	if (cellsMappingList.size()!=field.size())
	  {
	    FatalErrorIn
	      (
	       "setPODFields::int main"
	       )
	      <<"The cells mapping list size is not equal to field size "<<nl
	      <<"Mapping list size "<<cellsMappingList.size()<<nl
	      <<"Field size "<<field.size()<<nl
	      << exit(FatalError);
	  };

	Info<<"    Cells mapping: OK"<<nl<<endl;
	List<DynamicList< scalar > > timeCoeffsList(modesNumber);
	char sep = '\t';

	readTimeCoefficients(timeCoeffsFile, timeCoeffsList, sep);

        //work on files
        wordRe modesMaskName (modesFileName);
        if (wordRe::isPattern(modesFileName))
          {
            modesMaskName.compile(wordRe::REGEXP);
          }
        else
          {
            modesMaskName.compile(wordRe::LITERAL);
          }

        Foam::SortableList<fileName> modeDirFileList
          (
           readDir(inputModesDir, fileName::FILE)
           );

        int modeIndex = 0;

        forAll(modeDirFileList, i)
          {
            if(
               ! modesMaskName.match( modeDirFileList[i].name() ) 
               ||
               modeIndex >= modesNumber
               )
              continue;
            
            std::ostringstream ostr;
            ostr <<fieldName<<"_mode_"<<modeIndex;
            std::string converted = ostr.str();                  
            word modeName(converted);

            Info<<nl<<"Source File: "<<modeDirFileList[i]<<endl;
	    
            Foam::IFstream inputStream(
                                       inputModesDir
                                       /
                                       modeDirFileList[i]);
            
            Istream& modeStreamTest (inputStream);
            
            fieldType modeField
              (
               IOobject
               (
                modeName,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
                ),
               mesh,
               dimless,
               zeroGradientFvPatchScalarField::typeName
               );
	    
            forAll(modeField, celli)
              {
                label cellIndex = cellsMappingList[celli];
                const Type &value = pTraits<Type>(modeStreamTest);
                modeField[cellIndex] = value;
              }
	    
            if(saveFlag == "true") 
              {
                Info<<"Saving mode "<<modeIndex << " into "
                    <<mesh.time().timeName()<<nl<<endl;
                modeField.write();
              }
            if (relativeFlag == "true")
              {
                int timePoint(goalTime);
                forAll(field, celli)
                  {
                    field[celli] += 
                      timeCoeffsList[modeIndex][timePoint]
                      *
                      modeField[celli];
                  }
                
                int cellIndex = cellsMappingList[0];
                Info<<"Reconstruction time coeff = "
                      <<timeCoeffsList[modeIndex][timePoint]<<nl
                      <<"Mode value = "
                      <<modeField[0]<<nl
                      <<"Reconstructed value at point = "
                      <<field[cellIndex]<<nl
                      <<"Max reconstructed val = "
                      <<max(field)<<endl;           
              }
            modeIndex++;
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
      const labelList& mappingCells_;
      const settingsStruct& settings_;
      const fileName& dir_;

    public:

      iNew(const fvMesh& mesh, 
	   const labelList& mappingCells,
	   const settingsStruct& settings, 
	   const fileName& dir)
        :
	mesh_(mesh),
	mappingCells_(mappingCells),
	settings_(settings),
	dir_(dir)
      {}

        autoPtr<reconstructCellField> operator()(Istream& fieldValues) const
        {
            word fieldType(fieldValues);

            if
            (
	     !(
	       reconstructCellFieldType<scalar>
	       (fieldType, mesh_, mappingCells_, fieldValues, settings_, dir_)
	       || reconstructCellFieldType<vector>
	       (fieldType, mesh_, mappingCells_, fieldValues, settings_, dir_)
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

labelList getCellsMappingList
(
 const settingsStruct& PODSettings,
 const fvMesh& mesh,
 const IOdictionary& settingsDict,
 const fileName& caseDir
)
{
    DynamicList< label > cellsMappingList;

    word mode = PODSettings.mode;
    word relativeFlag = PODSettings.relativeFlag;   
    int modesNumber = PODSettings.modesNumber;
    word timeFileName = PODSettings.timeFileName;
    word modesFileName = PODSettings.modesFileName;

    wordRe modesMaskName (modesFileName);

    if (wordRe::isPattern(modesFileName))
      {
	modesMaskName.compile(wordRe::REGEXP);
      }
    else
      {
	modesMaskName.compile(wordRe::LITERAL);
      }

    Info<< "Setting field values using Proper Orthogonal Decomposition data" << endl;
    
    //Set up file names
    fileName inputModesDir
      (
       caseDir/"postProcessing"/"POD"/"modes"
       );

    fileName inputMappingDir
      (
       caseDir/"postProcessing"/"POD"/"mapping"
       );

    fileName timeMappingFile
      (
       caseDir/"postProcessing"/"POD"/"mapping"/"timeMapping"
       );

    //
    //Construction of the mapping list for cells (in parallel mode)
    // Pre-checks
    if ( modesNumber == 0
         || !isDir(inputModesDir)
         || ( relativeFlag == "false" && !isFile(timeMappingFile) ) 
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
      else
        {
          Info()<<nl<<"Error!"<<endl;
        }
      FatalErrorIn
	(
	 "setPODFields::int main()"
	 )
	<<"Invalid reconstruction options"<<nl
	<< exit(FatalError);
    }
    // Read directory entries into a list. Order makes sense

    Foam::SortableList<fileName> modeEntries
      (
       readDir(inputModesDir, fileName::FILE)
       );

    Info<<"    Found "
	<<modeEntries.size()
	<<" files in modes directory"<<endl;

    if(mode == "parallel")
      {
	// Read mapping entries
	// The mapping sequence should save order of processors
	Foam::SortableList<fileName> mapEntries
	  ( 
	   readDir(inputMappingDir, fileName::FILE) 
	    );
              
	Info<<"    Proceed in parallel mode. "
	    <<"Found "<<mapEntries.size()
	    <<" files in mapping directory"
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
	  <<"Unknown mode for reconstruction. "
	  <<"Use serial or parallel mode "<<nl
	  << exit(FatalError);
      }

    return cellsMappingList;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedMesh.H"

    Info<< "Reading setPODFieldsDict\n" << endl;

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

    fileName caseDir = runTime.path();
    settingsStruct globalPODSettings;

    readPODSettings(globalPODSettings, setPODFieldsDict);

    // Create mapping list
    labelList cellsMappingList = getCellsMappingList
      (globalPODSettings, mesh, setPODFieldsDict, caseDir);

    //now we got cellMappingList 

    if (setPODFieldsDict.found("reconstructFields"))
      {
        Info<< "Reconstructing field values" << endl;
        
	PtrList<reconstructCellField> reconstructedFieldValues
	  (
	   setPODFieldsDict.lookup("reconstructFields"),
	   reconstructCellField::iNew
	   (
	    mesh, 
	    cellsMappingList, 
	    globalPODSettings, 
	    caseDir
	    )
	   );
      }
    
    Info<< "\nEnd\n" << endl;

    return 0;
}

// ************************************************************************* //
