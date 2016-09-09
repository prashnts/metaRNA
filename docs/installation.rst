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

    git clone https://github.com/prashnts/metaRNA.git
    cd metaRNA
    python setup.py install

You can run the tests using the test-runner::

    python setup.py test

Browse the source code online at https://github.com/prashnts/metaRNA

Pre-requisite
-------------

`ViennaRNA <https://www.tbi.univie.ac.at/RNA/>`_ is required to compile
metaRNA C extensions. It is recommended to install ViennaRNA from source.
On Unix-like systems, it usually involves:

.. code-block:: bash

    wget 'http://www.tbi.univie.ac.at/RNA/download/sourcecode/2_2_x/ViennaRNA-2.2.10.tar.gz' -O viennarna.tar.gz
    mkdir viennarna
    tar -zxvf viennarna.tar.gz -C viennarna --strip-components=1
    cd viennarna
    ./configure
    make
    sudo make install

Usual build essentials (``automake``, ``autoconf``, ``gcc``) are required.

You can download the above package from `this link <https://noop.pw/etc/vienna224.tar.gz>`_
if the above link isn't accessible. To download and verify the SHA checksum:

.. code-block:: bash

    wget "https://noop.pw/etc/vienna224.tar.gz" -O vienna224.tar.gz
    echo "71a4c4704228fd01eb6e39415400a904d5240cef  vienna224.tar.gz" | shasum -c


Windows System
--------------

metaRNA hasn't been tested or built on Windows systems yet. Contributions
are welcome.

