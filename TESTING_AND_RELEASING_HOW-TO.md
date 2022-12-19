# Testing and releasing how-to

## Rationale
This markdown file describe the modalities of FijiRelax automatic tests and the release deployment pipeline. 

## Testing

### Automatic tests Rationale
FijiRelax tests are intended to support continuous integration of features and helping contributors to assess stability of the software through passes of code editing / bug correction / addition of new features, ....

### Dataset used for tests
Test execute on both synthetic data generated on the fly, and data stored in a repetition dataset directory, in src/test/resources/data/test_1/ . These data that can be used for a majority of FijiRelax functions.
See examples of accessing / generate these data in io.github.rocsg.fijirelax.test.Test_FijiRelaxPackage.java. 

### How and when to run automatic tests
Tests should be run before each commit, and have to be run before any deployment operation. That way, the building system of maven prevent to deploy any unstable version of FijiRelax.
To run automatic tests, you can run a maven build sequence that include testing. This action run a global lookup of classes/functions running tests and execute before deploying a local version of the jar file.
Examples of automatic tests can be found in io.github.rocsg.fijirelax.test.Test_FijiRelaxPackage.java. 

### How and when to add tests
Adding new features imply adding vulnerabilities. When contributing with new feature, please consider a test-driven development style, by first writing the tests, then writing the funcitons, and including the test of your new feature to the automatic JUnit test, by using a function decorator, or by adding the call to your function as a test_something function in io.github.rocsg.fijirelax.test.Test_FijiRelaxPackage.java
See examples of automatic tests in io.github.rocsg.fijirelax.test.Test_FijiRelaxPackage.java. 


## Deploying

### Compulsory deployment targets
Deployment imply to update these targets :
* Github repository : FijiRelax is open-source, thus the code of any release should be accessible to the world without any restriction. The release should be indicated in the github repository, by commiting a tag and then using it to create a release number.
* FijiRelax jar file, on ImageJ repository : FijiRelax jar is distributed through the ImageJ binaries repository "Fijiyama". People can use the ImageJ gui to check this repo to allow automatic installation ; hence they should access to the latest released build at any time.
* FijiRelax jar file, on Maven repository : FijiRelax is integrated as a maven artifact AKA io.github.rocsg.fijirelax ; similarly, the latest released version should be updated. Ask for romainfernandez06@gmail.com for committing to maven
* FijiRelax API, on javadoc : FijiRelax API is deployed on javadoc online, and possible contributors / extenders have to be able to access the latest update of the API. This can be done by logging on javadoc, and "asking" for a specific version of the FijiRelax API. That way, javadoc run a lookup of maven repository, and generate the doc on the fly, before storing it for being accessed by anyone.

### Suggested maintenance of version indicators at deployment time
Ideally, the "visual indicators" of version update that are non-automatically handled have to be edited to keep coherence.
* In io.github.rocsg.fijirelax.gui.FijiRelax_Gui.java, there are two String describing the release, which appear at execution in the GUI of FijiRelax : versionName="Handsome honeysuckle"; and timeVersionFlag="  Release time : 2022-10-24 - v4.0.3";
* In the README.md of the Github repo, there are badges linking to the last tag, the latest jar on maven, and the latest API	
* In the README.md of the ImageJ plugin page (https://imagej.net/plugins/fijirelax), there are the same badges



## Possible evolutions
As FijiRelax is part of ImageJ, it can be tough to design a github action, but a welcome evolution could be such a thing.


