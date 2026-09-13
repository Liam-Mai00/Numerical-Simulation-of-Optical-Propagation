#set text(
  font: "New Computer Modern",
  size: 12pt,
)

#set document(
  title: [Numpy Notes],
  author: "Liam Mai",
  date: auto,
  description: [Numpy notes],
)

#set heading(
  numbering: "1.1."
)

#figure(
  image("assets_notes/NumPy_logo.png")
)

#align(center)[
  #text(size:20pt)[
    #title()
    ]
    ]

#v(5em)
#align(center)[Liam Mai]
#pagebreak()
= Numpy Basics
Numpy is built around a homogeneous multidimensional array. _Homogeneous_ means that data within the array are of the same type e.g. int, double etc. The array is indexed using a tuple of non-negative integers. In addition, NumPy dimensions are called *axes*.

As an example, a 3D point in space could be ```python [1,2,1] ```it is 1 axis with 3 elements. However, with the next example:
```python [[1., 0., 0.],[0., 1., 2.]]```
The array consists of 2 axes (dimensions). The first axis has a length of 2 and the second with length 3.

Compared to MATLAB, Python uses row major indexing where, for 3D arrays, would be ```python (depth, row, column)``` instead of ```MATLAB (row, column, depth)```. Indexing arrays follow the same structure, so accessing the first and second layer of a matrix would be: ```python x[0,:,:]``` and ```python x[1,:,:]``` respectively.

