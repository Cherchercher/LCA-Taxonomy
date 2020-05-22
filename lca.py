# -*- coding: utf-8 -*-
"""LCA Taxonomy

Find the least common ancestor given two specie names
"""

import csv
import pandas as pd
from collections import defaultdict
import math
from io import BytesIO
from tarfile import TarFile
from urllib.request import urlopen

class TaxTree:
    def __init__(self, numNodes, level):
        self.numNodes = numNodes
        self.level = level
        self.tree = defaultdict(list)
        self.depth = [0]*(numNodes+1)
        self.parent = [[-1]*self.level for i in range(numNodes+1)]

    def dfs(self, cur, prev):
        self.depth[cur] = self.depth[prev] + 1;
        self.parent[cur][0] = prev
        for i in range(len(self.tree[cur])):
            if (self.tree[cur][i] != prev): 
                self.dfs(self.tree[cur][i], cur)
    
    def precomputeSparseMatrix(self):
        for i in range(self.level):
            for j in range(self.numNodes + 1):
                if (self.parent[j][i - 1] != -1): 
                    self.parent[j][i] = self.parent[self.parent[j][i - 1]][i-1] 

    def lca(self, u, v):
        if u == -1 or v == -1:
            return 'Not Found'
        if (self.depth[v] < self.depth[u]):
            u = u + v 
            v = u - v 
            u = u - v
        diff = self.depth[v] - self.depth[u]
        for i in range(self.level):
            if ((diff >> i) & 1) == 1:
                v = self.parent[v][i];
        if u == v:
            return u
        for i in range(self.level-1,-1,-1):
            if self.parent[u][i] != self.parent[v][i]:
                u = self.parent[u][i]
                v = self.parent[v][i]
  
        return self.parent[u][0]

    def addEdge(self, u, v):
        self.tree[u].append(v) 
        self.tree[v].append(u)

def find_name(idx, df):
    """
    lookup the tax_id in names.dmp (O(1))
    return tax_index if found, -1 if not found
    """
    name = df.loc[df['tax_id'] == idx,'name_txt']
    return -1 if len(name) == 0 else name.iloc[0]
    

def find_index(name, df):
    """
    lookup the tax_id in names.dmp (O(1))
    return tax_index if found, -1 if not found
    """
    idx = df.loc[df['name_txt'] == name,'tax_id']
    return -1 if len(idx) == 0 else idx.iloc[0]
    
def test_tree():
    tree = TaxTree(8, 1000)
    tree.addEdge(1,2) 
    tree.addEdge(1,3) 
    tree.addEdge(2,4) 
    tree.addEdge(2,5) 
    tree.addEdge(2,6) 
    tree.addEdge(3,7) 
    tree.addEdge(3,8)
    tree.dfs(1, 0)
    tree.precomputeSparseMatrix() 
    assert tree.lca(4, 7) == 1
    assert tree.lca(4, 6) == 2

def test(names):
    test_tree()
    assert find_index('root', names) == 1
    assert find_index('hi', names) == -1

def construct_tree_from_df(df, tree):
    for row in df.itertuples(index=True, name='Pandas'):
        tree.addEdge(row.parent_tax_id, row.tax_id)
    return tree

def construct_oneway_tree_from_df(df, tree):
    for row in df.itertuples(index=True, name='Pandas'):
        tree[row.parent_tax_id].append(row.tax_id)
    return tree

def read_files():
    print('reading files....')
    # or: requests.get(url).content
    names = pd.read_csv('names.dmp',sep='\t\|\t', lineterminator='\t\|\n')
    nodes = pd.read_csv('nodes.dmp',sep='\t\|\t', lineterminator='\t\|\n')
    #read df and get relevant cols from names
    names = names.iloc[:,0:2]
    names.columns = ['tax_id', 'name_txt']
    #read df and get relevant cols from nodes
    nodes = nodes.iloc[:,0:2]
    nodes.columns = ['tax_id', 'parent_tax_id']
    return names, nodes

def max_depth(tree, root):
    """
    :type root: Node
    :rtype: int
    """
    if len(tree[root]) == 0:
      return 0
    max_height = 0
    for child in tree[root]:
        max_height = max(max_height, max_depth(tree, child))
    return max_height+1

if __name__ == "__main__":
    # execute only if run as a script
    main()

def main():
    if len(sys.argv) != 2:
        print('please input two names seperated by space ')
        sys.exit()
    names, nodes = read_files()
    print('constructing tree....')
    num_nodes = names.shape[0]
    d_tree = defaultdict(list)
    one_way_tree = construct_oneway_tree_from_df(nodes, d_tree)
    max_level = max_depth(one_way_tree, 1) + 1
    tree = TaxTree(num_nodes, max_level)
    print('start constructing tree....')
    chunk_size = int(nodes.shape[0]/10)
    for start in range(0, nodes.shape[0], chunk_size):
        print("processing rows {} to {}".format(start, start+chunk_size))
        nodes_subset = nodes.iloc[start:start + chunk_size]
        tree = construct_tree_from_df(nodes_subset, tree)
    tree.dfs(1, 0)
    tree.precomputeSparseMatrix()
    print('finish constructing tree....')
    print('The lease common ancestor of {} and {} is {}'.format(sys.argv[0], sys.argv[1], find_name(tree.lca(int(find_index(sys.argv[0], names)), int(find_index(sys.argv[1], names))), names)))

#data vis
#new_df = pd.merge(names, nodes, on='tax_id', how='left')
#new_df[new_df.isna().any(axis=1)]
#data = []
#def construct_circular_tree_data_from_df(df):
#  for row in df.itertuples(index=True, name='Pandas'):
#    if row.name_txt == 'root':
#      data.append({"id":row.tax_id, "name": row.name_txt})
#    else:
#      data.append({"id":row.tax_id, "name": row.name_txt, "parent":row.parent_tax_id })
#construct_circular_tree_data_from_df(new_df)
#
#!pip install vega
#!jupyter nbextension install --sys-prefix --py vega  # not needed in notebook >= 5.3
#
#from vega import Vega
#def configure_vega_browser_state():
#  import IPython
#  display(IPython.core.display.HTML('''
#        <script src="/static/components/requirejs/require.js"></script>
#        <script>
#          window.outputs = [];
#          requirejs.config({
#            paths: {
#              base: '/static/base',
#              jquery: '//ajax.googleapis.com/ajax/libs/jquery/2.0.0/jquery.min',
#            },
#          });
#        </script>
#        '''))
#configure_vega_browser_state()
#Vega({
#  "$schema": "https://vega.github.io/schema/vega/v5.json",
#  "description": "Taxonamy Circular Tree Visualization.",
#  "width": 720,
#  "height": 720,
#  "padding": 5,
#  "autosize": "none",
#
#  "signals": [
#    {
#      "name": "labels", "value": True,
#      "bind": {"input": "checkbox"}
#    },
#    {
#      "name": "radius", "value": 280,
#      "bind": {"input": "range", "min": 20, "max": 600}
#    },
#    {
#      "name": "extent", "value": 360,
#      "bind": {"input": "range", "min": 0, "max": 360, "step": 1}
#    },
#    {
#      "name": "rotate", "value": 0,
#      "bind": {"input": "range", "min": 0, "max": 360, "step": 1}
#    },
#    {
#      "name": "layout", "value": "tidy",
#      "bind": {"input": "radio", "options": ["tidy", "cluster"]}
#    },
#    {
#      "name": "links", "value": "line",
#      "bind": {
#        "input": "select",
#        "options": ["line", "curve", "diagonal", "orthogonal"]
#      }
#    },
#    { "name": "originX", "update": "width / 2" },
#    { "name": "originY", "update": "height / 2" }
#  ],
#
#  "data": [
#    {
#      "name": "tree",
#      "values": data,
#      "transform": [
#        {
#          "type": "stratify",
#          "key": "id",
#          "parentKey": "parent"
#        },
#        {
#          "type": "tree",
#          "method": {"signal": "layout"},
#          "size": [1, {"signal": "radius"}],
#          "as": ["alpha", "radius", "depth", "children"]
#        },
#        {
#          "type": "formula",
#          "expr": "(rotate + extent * datum.alpha + 270) % 360",
#          "as":   "angle"
#        },
#        {
#          "type": "formula",
#          "expr": "PI * datum.angle / 180",
#          "as":   "radians"
#        },
#        {
#          "type": "formula",
#          "expr": "inrange(datum.angle, [90, 270])",
#          "as":   "leftside"
#        },
#        {
#          "type": "formula",
#          "expr": "originX + datum.radius * cos(datum.radians)",
#          "as":   "x"
#        },
#        {
#          "type": "formula",
#          "expr": "originY + datum.radius * sin(datum.radians)",
#          "as":   "y"
#        }
#      ]
#    },
#    {
#      "name": "links",
#      "source": "tree",
#      "transform": [
#        { "type": "treelinks" },
#        {
#          "type": "linkpath",
#          "shape": {"signal": "links"}, "orient": "radial",
#          "sourceX": "source.radians", "sourceY": "source.radius",
#          "targetX": "target.radians", "targetY": "target.radius"
#        }
#      ]
#    }
#  ],
#
#  "scales": [
#    {
#      "name": "color",
#      "type": "linear",
#      "range": {"scheme": "magma"},
#      "domain": {"data": "tree", "field": "depth"},
#      "zero": True
#    }
#  ],
#
#  "marks": [
#    {
#      "type": "path",
#      "from": {"data": "links"},
#      "encode": {
#        "update": {
#          "x": {"signal": "originX"},
#          "y": {"signal": "originY"},
#          "path": {"field": "path"},
#          "stroke": {"value": "#ccc"}
#        }
#      }
#    },
#    {
#      "type": "symbol",
#      "from": {"data": "tree"},
#      "encode": {
#        "enter": {
#          "size": {"value": 100},
#          "stroke": {"value": "#fff"}
#        },
#        "update": {
#          "x": {"field": "x"},
#          "y": {"field": "y"},
#          "fill": {"scale": "color", "field": "depth"}
#        }
#      }
#    },
#    {
#      "type": "text",
#      "from": {"data": "tree"},
#      "encode": {
#        "enter": {
#          "text": {"field": "name"},
#          "fontSize": {"value": 9},
#          "baseline": {"value": "middle"}
#        },
#        "update": {
#          "x": {"field": "x"},
#          "y": {"field": "y"},
#          "dx": {"signal": "(datum.leftside ? -1 : 1) * 6"},
#          "angle": {"signal": "datum.leftside ? datum.angle - 180 : datum.angle"},
#          "align": {"signal": "datum.leftside ? 'right' : 'left'"},
#          "opacity": {"signal": "labels ? 1 : 0"}
#        }
#      }
#    }
#  ]
#}
#)
