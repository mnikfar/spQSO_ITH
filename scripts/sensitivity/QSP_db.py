# -*- coding: utf-8 -*-
"""
Created on Tue May 26 20:11:19 2020
Save the result of a set of simulations to a database file
@author: Chang
"""

import sqlite3
from sqlite3 import Error

TABLE_NAME_PARAM = 'table_param' # 
TABLE_NAME_SIM = 'table_simulation' # 
TABLE_NAME_QSP = 'table_qsp'
TABLE_NAME_ABM = 'table_abm'


TYPE_STR = 'TEXT'
TYPE_FLOAT = 'REAL'
TYPE_INT = 'INTEGER'

class sim_result_db():
    # constructor
    def __init__(self, db_file):
        self.conn = None
        self.c = None
        self._connect(db_file)
        return
    # connect to database file; create one if not exist
    def _connect(self, db_file):
        try:
            self.conn = sqlite3.connect(db_file)
        except Error as e:
            print(e)
        self.c = self.conn.cursor()
        return
    # commit edit
    def commit(self):
        self.conn.commit()
        return
    # close connection
    def close(self):
        self.commit()
        self.conn.close()
        return
    # col_name and col_type: list of same length
    def create_table(self, table_name, col_name, col_type):
        cmd = 'CREATE TABLE IF NOT EXISTS {} (\n'.format(table_name)
        cmd += 'UID INTEGER PRIMARY KEY AUTOINCREMENT,\n'
        assert len(col_name) == len(col_type)
        cmd += ',\n'.join(['[{}] {}'.format(c, col_type[i]) for i, c in enumerate(col_name)])
        cmd += ')'
        self.c.execute(cmd)
        return
    # table info: cid, name, type, notnull, default, primary_key
    def get_table_info(self, table_name):
        return self.c.execute('PRAGMA TABLE_INFO({})'.format(table_name)).fetchall()
    # column names
    def get_colnames(self, table_name):
        return [c[1] for c in self.get_table_info(table_name)]
    # add one entry
    def add_entry(self, table_name, cols, values):
        col_name = ", ".join(['[' + c + ']' for c in cols])
        val_str = ", ".join(values)
        self.c.execute("INSERT OR IGNORE INTO {tn} ({col}) VALUES ({val})"\
                  .format(tn=table_name, col=col_name, val=val_str))
        return self.c.lastrowid
    # fetch results. colnames: list of column names (even if only one column); condition: string
    def fetch(self, table_name, col_names='*', condition=''):
        cols = '*'
        if col_names != '*':
            col_name_fix = ['[' + c + ']' for c in col_names]
            cols = ','.join(col_name_fix)
        cmd = "SELECT {} FROM {} {}".format(cols, table_name, condition)
        #print(cmd)
        self.c.execute(cmd)
        res = self.c.fetchall()
        return res
