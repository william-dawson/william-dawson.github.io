---
layout: post
title:  "What should be the output format for scientific codes?"
date:   2017-12-13 11:00:00 +0900
categories: blog best-practices
draft: true
---

Recently, I've been wondering what is the best output format
to use for a scientific code. By output, I don't mean just the raw data
that you have computed, but instead the real time output to the terminal that
happens as you perform the calculation.

For this questions, I came up with three features I would like to see in
a program output format.

1. Readable (that is, by a human)
2. Parsable (that is, by a computer)
3. Verifiable

Of these three, probably `Readable` is probably the most important. If we
didn't care about humans reading the output, we might as well skip writing
to the terminal in the first place. But users want to see that the program is
running as expected, and track progress on the calculation. `Readable` ensures
they can do this easily. Most programs you run will focus on having `Readable`
output.

But the program developer can't always anticipate all the data that a user
might want from a program. That's why `Parsable` is important. It should be
easily for the user to extract any piece of data about how a calculation was
performed and what its results were from the output file.

The least obvious point is `Verifiable`. Using just the output of a calculation,
we would like to be able automatically verify that a calculation was set up
correctly, and that it completed the calculation as expected.

Now let's consider some possible output formats and see how they match with my
criteria. As an example program, we'll consider the Hotelling's method program
from a previous blog post.

{% highlight python %}
from scipy.sparse.linalg import eigsh,norm
from scipy.sparse import identity
from scipy.io import mmread
from sys import argv

file_name = argv[1]
matrix = mmread(file_name)
dimension = matrix.shape[0]

largest_eigenvalue = eigsh(matrix, k=1, which='LM',
  return_eigenvectors=False)

inverse_mat = matrix * 1.0/(largest_eigenvalue[0]**2)

for i in range(0,100):
  inverse_mat = 2*inverse_mat - \
    inverse_mat.dot(matrix.dot(inverse_mat))
  norm_value = norm(inverse_mat.dot(matrix) - identity(dimension))
  if norm_value < 1e-8:
    break

{% endhighlight %}

## Spaghetti Output
Probably the most common output format is... well no format at all. Let's
imagine what that might look like with our program.

```
############################################################################
/  |      /         / / /                     /|/|      /    /            |
(___| ___ (___  ___ ( (    ___  ___  ___      ( / | ___ (___ (___  ___  ___|
|   )|   )|    |___)| | | |   )|   )|___      |   )|___)|    |   )|   )|   )
|  / |__/ |__  |__  | | | |  / |__/  __/      |  / |__  |__  |  / |__/ |__/
############################################################################
Reading input...
Output read! Input file=S_matrix.mtx
Computing initial guess matrix...
Largest eigenvalue:1.86124928636
Initial guess computed!
Running algorithm
1.57802206403
1.2555463121
0.966656384587
0.652685328614
0.312962675887
0.0782546324211
0.00561087257415
3.11928836745e-05
Program finished, with convergence of: 9.72911618378e-10
Thank you for using our program
```

Ugly, but at least it's friendly. And of course I couldn't forget the giant
ascii header that we've all seen before. This output gives us enough output to
know what stage we are at in the calculation, so a user can easily track the
progress of the code. We could try to be even more pretty with our
output. How about this:


```
# My Hotelling's Method Code
Compute a matrix inverse using hotelling's method

## Reading input...
Output read! Input file=S_matrix.mtx

## Computing initial guess matrix...
Largest eigenvalue: 1.86124928636
Initial guess computed!

## Running algorithm
Iterations
0. 1.57802206403
1. 1.2555463121
2. 0.966656384587
3. 0.652685328614
4. 0.312962675887
5. 0.0782546324211
6. 0.00561087257415
7. 3.11928836745e-05

## Summary
Program finished, with convergence of: 9.72911618378e-10
Thank you for using our program
```

That's right, it's markdown. If the text alone isn't easy enough to read, you
could easily convert it into an html website even. Maybe I could just run my
programs instead of writing blog posts.

So what are the downsides of this approach? Well mainly it is difficult to parse.
Not impossible, but you basically have to write a custom parser that knows
where in the file all the values you want are. And changes to code could
competely break your parser.

## Formatted Output

Let's consider some alternative formats. First, what about XML:

```
<?xml version="1.0" encoding="ISO-8859-1" ?>
<Hotelling_Method>
  <message>Compute a matrix inverse using Hotelling's method</message>
  <reading_input>
    <input_file>S_matrix.mtx</input_file>
  </reading_input>
  <initial_guess>
    <largest_eigenvalue>1.86124928636</largest_eigenvalue>
  </initial_guess>
  <run_algorithm>
    <value>1.57802206403</value>
    <value>1.2555463121</value>
    <value>0.966656384587</value>
    <value>0.652685328614</value>
    <value>0.312962675887</value>
    <value>0.0782546324211</value>
    <value>0.00561087257415</value>
    <value>3.11928836745e-05</value>
  </run_algorithm>
  <summary>
     <convergence>9.72911618378e-10</convergence>
     <message>Thank you for using our program</message>
  </summary>
</Hotelling_Method>
```

The Qbox[1] electronic structure package, for example, uses XMl.
One nice thing about using XML is that there are many libraries that could be
used to parse it. Additionally, the tree like organization of XML means
a parse can quickly get to the values of interest. But as a downside, this
XML can be a little difficult to read for a human.

Currently, I am working on collaborations with the developers of BigDFT[2].
BigDFT's choice of output is YAML:

```
Hotelling_Method:
  message:"Compute a matrix inverse using Hotelling's method"
  input:
    input_file:S_matrix.mtx
  initial_guess:
    largest_eigenvalue:1.86124928636
  run_algorithm:
    0: 1.57802206403
    1: 1.2555463121
    2: 0.966656384587
    3: 0.652685328614
    4: 0.312962675887
    5: 0.0782546324211
    6: 0.00561087257415
    7: 3.11928836745e-05
  summary:
     convergence:9.72911618378e-10
     message:"Thank you for using our program"
```

In my opinion, YAML is less verbose, which makes it much easier for the user to
read. And like XML, there are many libraries for parsing YAML automatically.
JSON would also make for a nice alternative.

## Verifiable

The one benefit I see for XML over YAML is the notion of XML schema. Being able
to automatically validate the output file can ensure the code remains parsable
despite the introduction of new features.
Furthermore, schema can help a user get a high level view of where output values
are stored.

I would like to go one step further with the concept of verifiable though. One
amazing feature about BigDFT, is that every output file is also an input file!
That means you should be able to truly verify output by rerunning the
calculation on the output file. This is a fantastic guard against data
tampering. You might hope to go even one step further and have the output file
contain all the information about how the program was compiled, and which
version of the code was used.

In the future, I hope to apply these principles to my program, and NTPoly of
course. Indeed, already it outputs in a YAML format. Unfortunately, I am still
struggling with how to make the output
printer best interact with the calling program. Perhaps my progress will make
for a good future blog post.

> [1] http://www.qboxcode.org/

> [2] http://bigdft.org/Wiki/index.php?title=BigDFT_website
