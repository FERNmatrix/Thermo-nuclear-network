Once unzipped and untarred, position yourself in the directory StochasticElements.  Compile the java with

```
javac *.java
```

and then execute 

```
java ElementMaker
```

This should launch a GUI to control the calculations.  There is a Help button on the GUI with information, but if you get the GUI to launch you can do a sample calculation by

1.  Click the Select Isotopes button (lower right)

2.  In the resulting popup, Save Changes to select the default network isotopes.

3.  In the resulting popup,  Save Changes to select the default reactions for the network.

4.  Some of the isotope blocks should now be purple (indicating that they have been selected for the network).  Shift-click on one of the purple blocks and in the popup click Save Changes to select the default initial abundances.

5.  Click the Parameters button on the right.  In the resulting popup, select Save Changes to select the default values of the parameters.

6.  Click the Integrate button on the right.

7.  This should generate a lot of output to the screen and at the end of the calculation launch a plot window with the results of the network integration.  I usually redirect the screen output to a file.  For example, in a unix/linux shell

```
java ElementMaker | tee temp.txt
```

will output to the screen and also to the file temp.txt.

There are many java files but the main file containing the algorithm is StochasticElements.java.
