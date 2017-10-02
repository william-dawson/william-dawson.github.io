---
layout: post
title:  "MPI Parallel File I/O For Text Files"
date:   2017-09-01 11:00:00 +0900
categories: blog method
draft: true
---

Today, I'm going to teach you how to perform parallel file I/O on text files.
Ideally, you shouldn't do this. Instead, the performance (and ease of coding)
is superior when performing I/O on binary files. Additionally, there are a
number of parallel I/O libraries[1][2] out there that can easily read and write
complex data.

However, there are some benefits for text files. A lot of the publicly available
data is in text formats. Text files are much easier for users to understand,
and are portable across different types of machines. You can also imagine
the benefits of writing data in standard text formats such as JSON or XML.

In those cases, we have to do parallel file I/O. I will show an example of this
here using the
[matrix market format](http://math.nist.gov/MatrixMarket/formats.html). We will
restrict ourselves to the simplest case: a dense matrix of real numbers. The
file format looks like this:

> Some comment lines
> matrix_rows matrix_columns
> data_11
> data_12
> etc

We can easily generate matrix market test files using python:
{% highlight python %}
from numpy.random import rand
from scipy.io import mmwrite
from sys import argv

matrix_size = int(argv[1])
output_file = argv[2]

matrix = rand(matrix_size, matrix_size)
mmwrite(output_file, matrix)
{% endhighlight %}

The main challenge you will face when doing MPI parallel file I/O on text files
is that for MPI you can't specify which lines to read. Instead, you can only
tell it which characters you want to read. This is a problem, because we don't
know how many characters per line there are before reading the file. And of
course it might be variable across lines. So each processor doesn't know which
characters to read in.

Here will be the main outline of our solution:

1. Read the header information.
2. Determine the size of the file.
3. Split the file size up evenly amongst processors.
4. Add some buffer size to read.
5. Read in parallel.
6. Pull out a subset of the read data.

The trick is in step 4 and 6. If we only read at character offsets, one
processor might find itself reading from the middle of a line, while another
processor finds itself ending in the middle of a line. So what we do is we
add some buffer, so that every processor reads at least one extra line. Then
if you're in the middle of a line, you skip that line. And if you end in the
middle of a line, you use that buffer to get the rest of the line.

Let's go over this in detail now. First, our basic setup:
{% highlight python %}
# Libraries
from sys import argv
from mpi4py import MPI
import array
import numpy

# MPI Setup
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
total_processors = comm.Get_size()

# Process Parameters
in_file = argv[1]
out_file = argv[2]
{% endhighlight %}

Next, we'll have the root processor parse the header information using a
sequential read.

{% highlight python %}
# Root Gets Basic Matrix Information
rows = 0
columns = 0
header_text = ""
if rank == 0:
  with open(in_file, 'r') as ifile:
    line = ifile.readline()
    header_text = line
    while line[0] == "%":
      line = ifile.readline()
      header_text = header_text+line
    split = line.split()
    rows = int(split[0])
    columns = int(split[1])

# Broadcast to everyone else
rows = comm.bcast(rows, root=0)
columns = comm.bcast(columns, root=0)
header_length = comm.bcast(len(header_text), root=0)
{% endhighlight %}

Now we open up the file, and get the total file size.

{% highlight python %}
# Parallel Read
fh = MPI.File.Open(comm, in_file, MPI.MODE_RDONLY)
total_file_size = MPI.File.Get_size(fh)
{% endhighlight %}

We now use the information about the total file size and the header to split
up the remaining parts of the file.

{% highlight python %}
## Compute Offsets and Amount of Data To Read
length_minus_header = total_file_size - header_length
max_line_size = 81
local_read_size = max(max_line_size, int(length_minus_header/total_processors))
local_read_size_plus_buffer = local_read_size + max_line_size
local_start_read = local_read_size*rank + header_length

## Handle small files
if local_start_read > total_file_size:
  local_start_read = 0
  local_read_size = 0
  local_read_size_plus_buffer = 0
elif local_start_read + local_read_size_plus_buffer > total_file_size:
  local_read_size = total_file_size - local_start_read
  local_read_size_plus_buffer = local_read_size
{% endhighlight %}

An even way to split up the remaining file is to take the total file size minus
the header size, and divide it by the number of processors. But remember we
also need some buffer space in case we finish reading halfway through a line.
So we add the maximum line size to our original read size. Finally, this portion
of code has to take care of the edge case where we have way too many processors
for a small file. In that case, we make sure the first few processors read at
least one line, and the remaining read nothing. You might also consider
switching to a serial read in this case.

Now we do the actual read. There is some extra, mpi4py specific code here that
comes from our need to set the datatype right for the read.

{% highlight python %}
local_buffer = numpy.empty(local_read_size_plus_buffer,dtype="c")
fh.Read_at_all(local_start_read, [local_buffer, MPI.CHAR])

local_data_string = b''.join(local_buffer).decode('utf-8')
fh.Close()
{% endhighlight %}

The last step is to get rid 

{% highlight python %}
# Remove Duplicate Data
start_char = 0
if rank > 0:
  while(local_data_string[start_char] != '\n'):
    start_char = start_char+1

end_char = local_read_size
while(end_char < local_read_size_plus_buffer and \
  local_data_string[end_char] != '\n'):
  end_char = end_char + 1

local_data_string = local_data_string[start_char:end_char]
{% endhighlight %}

> https://support.hdfgroup.org/HDF5/
> https://www.unidata.ucar.edu/software/netcdf/
