                         VecVis

Summary
-------
VecVis is a special purpose visualisation tool designed
for the effective visualisation of 2 dimensional, time
dependent vector fields with complex dynamics. VecVis 
utilises the Image-Base Flow Visualisation (IBFV) 
originally developed by Jarke van Wijk (2002).

Authors
-------

David J. Warne [1,2]
    Email: david.warne@qut.edu.au

Joseph Young [1]
    Email: j.young@qut.edu.au

[1] High Performance Computing and Research Support,
    Information Technology Services Department,
    Technology,Information, and Learning Support Division,
    The Queensland University of Technology,
    Brisbane QLD 4001, Australia

[2] School of Electrical Engineering and Computer Science,
    Faculty of Science and Engineering,
    The Queensland University of Technology,
	Brisbane QLD 4001, Australia


The Latest Version
------------------
The latest version of VecVis, including source code can 
be found on GitHub at https://

Documentation
-------------
This program is very much research code, thus detailed 
documentation is not currently practical. Please refer to
the Example section for help, or conatct the lead author 
at david.warne@qut.edu.au. 

Compiling the Code
------------------
If the Pre-built binaries do not work for your flavour of 
Operating system (Or you would prefer to compile it for
performance reasons) then please follow the compilation
steps for your operation system

Linux: If you are using Linux, you probably don't need 
       the following steps, but just in case :)...

	1. Ensure that freeglut is installed, if not then 
	   install it. For example on Fedora,

        $ yum install freeglut
    
    2. Run make

        $ make

    3. Thats it!

Windows: For now this has only been tested with MS Visual
         Studio 2010.
	1. Start-up MS Visual Studio 2010.

	2. Select File->Open Project and select the VecVis.sln file.

	3. Click ``Build''.

	4. Done!

Mac:

    1. Sadly, this software has not been tested for Mac.
       However, in theory it should be similar to Linux.
       Some include and library paths may need to be 
       modified, but the code should be portable... In
       theory.

``In theory, there is no difference between theory and
practice. But, in practice, there is'' ~ Albert Einstein

Licensing
---------
Please see the LICENSE.txt file.

Example
-------

Change to the bin/ directory

$ cd bin/

To view commandline options use the -h option,

$ VecVis -h

To check the version use the -v option,

$ VecVis -v

Run the Example script,

$ RunExample.sh

or the *.bat file on a Windows machine

Once VecVis is running you cen use the ? key to view the
list of keystroke commands. Dye and observations points 
can be selected with a mouse left click and moved with a
mouse right click.

The Example script uses the data set located in SampleData
which is a turbulent flow dataset built from numerical 
simulations. There are two files the *.ans and the *.flw. 
The *.ans is an annotations file, it holds text labels and
colourmap data (for scalar data overlays). The *.flw file
contains the flow/scalar data a polygonal mesh. It is a 
simple ASCII text format using a # to specify data types.
