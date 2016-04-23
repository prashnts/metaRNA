.. _installation:

Installing
==========

Using PyPi
----------

metaRNA can be installed very easily using `pip <http://www.pip-installer.org/en/latest/index.html>`_.

.. code-block:: bash

    pip install metarna



Using git
---------

If you want to run the very latest, feel free to pull down the repo from github
and install by hand.

.. code-block:: bash

    git clone https://github.com/PrashntS/metaRNA.git
    cd metaRNA
    python setup.py install

You can run the tests using the test-runner::

    python setup.py test

Browse the source code online at https://github.com/PrashntS/metaRNA

Pre-requisite
-------------

`ViennaRNA <https://www.tbi.univie.ac.at/RNA/>`_ is required to compile
metaRNA C extensions. It is recommended to install ViennaRNA from source.
On Unix-like systems, it usually involves:

.. code-block:: bash

    wget "http://www.tbi.univie.ac.at/RNA/download/package=viennarna-src-tbi&flavor=sourcecode&dist=2_2_x&arch=src&version=2.2.4" -O viennaRNA.tar.gz
    tar -zxvf ViennaRNA-2.2.5.tar.gz
    cd ViennaRNA-2.2.5
    ./configure
    make
    sudo make install

Usual build essentials (``automake``, ``autoconf``, ``gcc``) are required.

Windows System
--------------

metaRNA hasn't been tested or built on Windows systems yet. Contributions
are welcome.

