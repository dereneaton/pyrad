pyrad
=====

Assembly and analysis of RADseq data sets


Tutorials
---------

Detailed information and a number of example tutorials are 
available `here <http://www.dereneaton.com/software/pyrad/>`_.    


Downloads
---------

Stable release versions can be downloaded `here <https://github.com/dereneaton/pyrad/releases>`_, or you can clone the current development version using git:

::

    git clone https://github.com/dereneaton/pyrad.git



Installation (As of v.3.0.6 and newer)
^^^^^^^^^^^^^^^^^
With the following you can install pyrad so that it is callable as an executable from anywhere on your machine. If you have pip, then the second option will also install the dependencies numpy and scipy:

::

    cd pyrad
    sudo pip install .
    pyrad -h

Or 

::

    cd pyrad
    sudo python setup.py install  
    pyrad -h

Or alternatively, without having to install you can simply call pyRAD.py from its location using python:

::
    
    python pyrad/pyrad/pyRAD.py -h


Python requirements
^^^^^^^^^^^^^^^^^^^
You will need the following two Python dependencies with `pyrad`.

 * numpy
 * scipy

In addition to the programs  

 * `muscle <http://www.drive5.com/muscle/downloads.htm>`_
 * `vsearch <https://github.com/torognes/vsearch>`_

Example usage (see tutorials)
^^^^^^^^^^^^^^^^^^^
::

    >>> pyrad -n  
        new params.txt file created


    >>> pyrad -p params.txt 



Licence
-------
GPLv3  


Authors
-------

`pyrad` was written by `Deren Eaton <deren.eaton@yale.edu>`_.
