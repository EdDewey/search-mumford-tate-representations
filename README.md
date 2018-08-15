# search-mumford-tate-representations

Suppose you have a polarized Hodge structure and want to know its Mumford-Tate
group.  It is much easier to find the Hodge numbers, and if your Hodge
structure came from topology you probably know whether it is odd or even.  

For small-rank examples this constrains the possible Mumford-Tate group quite a
bit.  The file searchRepresentations.py contains  code to ennumerate the
possible Cartan types of the Mumford-Tate group, given only the parity of the
bilinear form and the Hodge numbers. 

The Jupyter notebook example.ipynb gives a usage example and sketches how the
enumeration works.  
