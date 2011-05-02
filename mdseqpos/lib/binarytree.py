#! /usr/bin/env python

"""A module for binary trees."""

import os

os.environ['DJANGO_SETTINGS_MODULE'] = 'mdseqpos.settings'

from django.template import loader, Context

class BinaryTree:
    """A class for binary trees."""
    
    ##################
    # initialization #
    ##################

    def __init__(self, v=None, l=None, r=None):
        """Initialize tree node."""
        self.value = v   # value for this tree node
        self.lchild = l  # pointer to the left child of this tree node
        self.rchild = r  # pointer to the right child of this tree node
        # check that left and right children are both tree nodes
        if self.lchild is not None and not isinstance(self.lchild, BinaryTree):
            raise TypeError, "left child %s is not of %s" % (self.lchild, BinaryTree)
        if self.rchild is not None and not isinstance(self.rchild, BinaryTree):
            raise TypeError, "right child %s is not of %s" % (self.rchild, BinaryTree)
        # set size of the tree
        self.set_size()
    
    def set_size(self):
        """Calculate the size of the tree."""
        self.size = 1
        if self.lchild is not None:
            self.size += self.lchild.size
        if self.rchild is not None:
            self.size += self.rchild.size
    
    ######################
    # consistency checks #
    ######################

    def check_type(self, t):
        """Check that values of all tree nodes are of type t."""
        if self.value is not None and not isinstance(self.value, t):
            raise TypeError, "value of %s is not of %s" % (self, t)
        if self.lchild is not None:
            self.lchild.check_type(t)
        if self.rchild is not None:
            self.rchild.check_type(t)
    
    ##################
    # output methods #
    ##################
    
    def map(self, function, kwargs={}):
        """Apply function to each node and return a tree of the results."""
        v = function(self.value, **kwargs)
        if self.lchild is None:
            l = None
        else:
            l = self.lchild.map(function, kwargs)
        if self.rchild is None:
            r = None
        else:
            r = self.rchild.map(function, kwargs)
        return BinaryTree(v, l, r)
    
    def nested_list(self):
        """Return the tree as a nested list."""
        children = []
        if self.lchild is not None:
            children.append(self.lchild.nested_list())
        if self.rchild is not None:
            children.append(self.rchild.nested_list())
        return [self.value, children]
    
    def to_html(self, binarytree_django_template='binarytree.html'):
        """Return the tree in HTML as a nested, unordered list."""
        def add_binarytree_connector(binarytree_node):
            """Add surrounding tags to each node for drawing lines that depict the structure of tree."""
            t = loader.get_template('binarytree-connector.html')
            c = Context({'binarytree_node': binarytree_node})
            return t.render(c)
        t = loader.get_template(binarytree_django_template)
        c = Context({'binarytree_as_nested_list': self.map(add_binarytree_connector).nested_list()})
        return t.render(c)

    def to_json(self, dst_dir="results/", img_dir="img/"):
        """Return the tree as json"""
        if self.value is None:
            return "{}"
        
        children = [];
        if self.lchild is not None:
            children.append(self.lchild.to_json(dst_dir, img_dir))
        if self.rchild is not None:
            children.append(self.rchild.to_json(dst_dir, img_dir))

        if len(children) > 0:
            if len(children) == 1:
                child_str = "'children':["+children[0]+"]"
            else: #there are two
                child_str = "'children':["+children[0]+","+children[1]+"]"
        else:
            child_str = "'children':[]"
            
        return "{'node': "+self.value.to_json(dst_dir, img_dir)+", "+child_str+"}"

    def to_new_html(self, django_template='mdseqpos_out.html', dst_dir='results/', img_dir='img/'):
        """Return the tree as HTML output--injecting the data as JSON with this
        structure {'node': new MotifHelper(...), 'children': [{NODES}]}
        """
        tree_to_json = self.to_json(dst_dir, img_dir)
        #print tree_to_json
        t = loader.get_template(django_template)
        c = Context({'tree_to_json': tree_to_json})
        return t.render(c)
    
    #####################
    # traversal methods #
    #####################
    
    def preorder(self):
        """Return a generator that traverses the tree in preorder."""
        if self.value is not None:
            yield self.value
        if self.lchild is not None:
            for value in self.lchild.traverse():
                yield value
        if self.rchild is not None:
            for value in self.rchild.traverse():
                yield value
    
    def inorder(self):
        """Return a generator that traverses the tree in inorder."""
        if self.lchild is not None:
            for value in self.lchild.traverse():
                yield value
        if self.value is not None:
            yield self.value
        if self.rchild is not None:
            for value in self.rchild.traverse():
                yield value
    
    def postorder(self):
        """Return a generator that traverses the tree in postorder."""
        if self.lchild is not None:
            for value in self.lchild.traverse():
                yield value
        if self.rchild is not None:
            for value in self.rchild.traverse():
                yield value
        if self.value is not None:
            yield self.value
    
    def traverse(self, order='pre'):
        """Return a generator that traverses the tree in the given order.
        
        order can take on the following values:
        'pre'  <-- traverse tree in preorder
        'in'   <-- traverse tree in inorder
        'post' <-- traverse tree in postorder"""
        if order == 'pre':
            return self.preorder()
        elif order == 'in':
            return self.inorder()
        elif order == 'post':
            return self.postorder()
        else:
            raise ValueError, "order must be 'pre', 'in', or 'post'"
