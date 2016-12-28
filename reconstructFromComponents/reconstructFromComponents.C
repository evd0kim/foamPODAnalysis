/*---------------------------------------------------------------------------* \
License
    This file is designed to work with OpenFOAM.

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

Application
    interpolateVolumeField

Description
    Interpolate vol*Field to surface*Field.

Author
    wyldckat@github

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "volFields.H"

#include "OSspecific.H"

#include "IOstreams.H"
#include "IStringStream.H"
       
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{

    timeSelector::addOptions();

    argList::addNote
    (
        "Reconstructing volumeVectorField from it's components"
    );
    argList::addOption
    (
        "vector",
        "volumeVectorField"
    );

#   include "addRegionOption.H"

#   include "setRootCase.H"
#   include "createTime.H"
    
    fileName vectorFieldName;

    bool vectorBool = args.optionReadIfPresent
      (
       "vector", vectorFieldName
       );

    if(!vectorBool)
      {
        FatalErrorIn(args.executable())
          << "In the current version the type of volume field should be -scalar or -vector."
          << exit(FatalError);
      }

#   include "createNamedMesh.H"

    Info<< "Time: " << runTime.timeName() << endl;

    mesh.readUpdate();

    if(vectorBool)
      {
        Info<< nl 
            << "Reading volume field " << vectorFieldName<<nl 
            << endl;
        
        volVectorField volVecField
          (
           IOobject
           (
            vectorFieldName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
            ),
           mesh
           );

        List<word> componentNames(IStringStream("(x y z)")());
        
        if (volVecField.size())
          {
            int idx = 0;
            for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
              {
                fileName compFieldName
                  (
                   vectorFieldName+componentNames[idx]
                   );
                
                Info<<"    Component "
                    <<cmpt<<", Projection "
                    <<componentNames[idx]<<nl
                    <<"    Field component "
                    <<compFieldName
                    <<endl;
                
                volScalarField volScaField
                  (
                   IOobject
                   (
                    compFieldName,
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                    ),
                   mesh
                   );
                
                volVecField.replace(cmpt, volScaField);
                idx++;
              }
          }
        volVecField.write();
      }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
