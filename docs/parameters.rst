.. _parameters:

Calculation Parameters
======================

miRanda algorithm used by metaRNA accepts the following optional parameters.

+-------------------------+----------+---------+-------------------------------------------------------------------+
| Parameter               | Type     | Default | Description                                                       |
+=========================+==========+=========+===================================================================+
| scale                   | float    | 4.0     | The 5' miRNA scaling parameter.                                   |
+-------------------------+----------+---------+-------------------------------------------------------------------+
| strict                  | int      | 0       | Perform a Strict Seed search when set to 1.                       |
+-------------------------+----------+---------+-------------------------------------------------------------------+
| gap_open                | float    | -9.0    | Gap-open Penalty                                                  |
+-------------------------+----------+---------+-------------------------------------------------------------------+
| gap_extend              | float    | -4.0    | Gap-extend Penalty                                                |
+-------------------------+----------+---------+-------------------------------------------------------------------+
| score_threshold         | float    | 50.0    | Score Threshold for reporting hits                                |
+-------------------------+----------+---------+-------------------------------------------------------------------+
| energy_threshold        | float    | 1.0     | Energy Threshold for reporting hits                               |
+-------------------------+----------+---------+-------------------------------------------------------------------+
| length_5p_for_weighting | int      | 8       | The 5' sequence length to be weighed except for the last residue. |
+-------------------------+----------+---------+-------------------------------------------------------------------+
| temperature             | int      | 30      | Used while calculating Free Energy                                |
+-------------------------+----------+---------+-------------------------------------------------------------------+
| alignment_len_threshold | int      | 8       | Minimum alignment.                                                |
+-------------------------+----------+---------+-------------------------------------------------------------------+


Passing Parameters
------------------

The parameters are passed as keyword arguments.

.. code:: python

    targets = scan(gene_sequence, mirna_sequence, scale=5.0, strict=1)

