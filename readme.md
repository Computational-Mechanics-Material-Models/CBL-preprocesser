# CBL Preprocessor

&#x20;[![BSD License](http://www.projectchrono.org/assets/logos/chrono-bsd.svg)](https://github.com/Cusatis-Computational-Services/CBL-preprocesser/blob/main/LICENSE)

This is a GitLab repository for the development and implementation of a preprocessor in [FreeCAD](https://www.freecadweb.org/) for the [Connector Beam Lattice Model](https://github.com/Cusatis-Computational-Services/CBL-chrono) in [Project Chrono](https://www.projectchrono.org) software. The repository and wiki are in development and may be incomplete or include mistakes.


## Desktop Installation and Setup

Follow the below instructions to get set up with the Wood Workbench for CBL. These instructions were written with Windows users in mind. It is also possible to install on Mac or Linux.

Please note the versions of each piece of software. Newer or alternate versions may work, but have not been tested and verified for compatibility.

<details>

<summary>Step 1: Install FreeCAD</summary>

Install the latest version of FreeCAD (use at least version 0.20.2). The download is available for free:

[https://www.freecadweb.org/downloads.php](https://www.freecadweb.org/downloads.php)

</details>

<details>

<summary>Step 2: Install Git Client</summary>

Any Git client can be used to push and pull from the GitHub. Two well-supported options are Github Desktop or Sourcetree.

[https://desktop.github.com/download/](https://desktop.github.com/download/)

[https://www.sourcetreeapp.com/](https://www.sourcetreeapp.com/)

Instructions below will cover these two options, though any client may work, as well as cloning via command line.

</details>

<details>

<summary>Step 3: Pull Wood Workbench</summary>

We recommend pulling the GitHub directly into the FreeCAD workbench directory. Otherwise if you pull to another location then you will need to copy the pulled files to the appropriate directory. On a Windows system this is likely under **C:\Users\\\<usr>\AppData\Roaming\FreeCAD\Mod**, where **\<usr>** is the system user.

1. Github Dekstop
* Open Github Desktop
* Select **File** > **Clone repository**
* Select "**URL**"
* Enter https://github.com/Cusatis-Computational-Services/CBL-preprocesser as the source URL.
* Paste the FreeCAD Mod path for the local path (or click **Choose** to browse to that path). If not already existing, a new folder called CBL-preprocessor should be created in this location.
* Click "**Clone**"

2. Sourcetree
* Open Sourcetree
* Select **File** > **Clone / New...**
* Select "**Remote**" and "**Add an account...**"
* For "Hosting Service" select "**GitHub**". For "Authentication" select "**OAuth**"
* Click on "**Refresh OAuth Token**" and login to GitHub and allow Sourcetree in the browser window that opens
* Click **Ok** in Sourcetree. Then "**CBL-preprocesser**" should populate on the right side of the window. If it doesn't, you may need to click refresh.
* Select "**CBL-preprocesser**" and click "**Clone**"
* When filling out the clone window, be sure to place the repository in the FreeCAD Mod path as noted above.
* Click "**Clone**"



</details>

<details>

<summary>Step 4: Supporting Packages</summary>

The Wood Workbench also requires the Triangle python package installed to the FreeCAD environment. 

* Navigate to **C:\Program Files\FreeCAD 0.21\bin** and open a terminal at this location.
* Run ``.\python.exe -m pip install triangle``

</details>

<details>

<summary>Step 5: Test Installation</summary>

Verify that everything is installed properly by opening FreeCAD. 

* Select the Wood Workbench from the workbench dropdown menu. 
* Select the particle icon (this should be the only workbench-specific button).
* Run the workbench in debug mode by scrolling to the bottom without changing any parameters and click **Generate Model**. 

If no erros appear during this process, the workbench is installed correctly.

</details>


## License

The CBL Preprocessor workbench application is licensed free and open source under the [BSD 3-Clause “New” or “Revised” License](https://choosealicense.com/licenses/bsd-3-clause/). One reason for Project Chrono's popularity is that its users are granted the freedom to modify and redistribute the software and have a right of continued free use, within the terms of the BSD license.


## Citing this Work

A proper citation for theses developments is still to come. In the meantime, please reference this GitHub in all publications.
