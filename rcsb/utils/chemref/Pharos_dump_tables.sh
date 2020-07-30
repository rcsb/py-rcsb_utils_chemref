#!/bin/bash
# File: Pharos_dump_compound_tables.sh
#
# Dump tables from the Pharos schema (TCRDv6.7) containing ChEMBL references...
#
mysql --user=$MYSQL_DB_USER --password=$MYSQL_DB_PW  -e "use tcrd6; select * from drug_activity;" >  drug_activity.tdd
mysql --user=$MYSQL_DB_USER --password=$MYSQL_DB_PW  -e "use tcrd6; select * from cmpd_activity;" >  cmpd_activity.tdd
