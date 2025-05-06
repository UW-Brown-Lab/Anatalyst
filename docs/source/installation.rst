Installation
============

Requirements
-----------

SC Pipeline requires:

TBD


Basic Installation with Docker
-----------------

You can use the provided Docker container which includes all dependencies:
*Right now this only contains the software dependencies, none of the actual pipeline code. WIP*

.. code-block:: bash

    docker pull ghcr.io/wblashka/downstream-scrnaseq:latest
    docker run -v /path/to/your/data:/data -it ghcr.io/wblashka/downstream-scrnaseq:latest