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
Using pip you can install pyrad so that it can be called from anywhere on your machine:

::

    cd pyrad
    sudo pip install .
    pyrad -h    

or alternatively, without installation you can simply use python to call pyrad from its location:

::

    python pyrad/pyRAD.py -h


Python requirements
^^^^^^^^^^^^^^^^^^^
These will be automatically installed alongside `pyrad` when using `pip`.

 * numpy
 * scipy


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
