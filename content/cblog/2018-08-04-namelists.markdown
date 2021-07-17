Title: A Program For Managing Namelists
Date: 2018-08-04 11:00:00 +0900
Category: Lessons

Recently, I like to stay in the office late once or twice a week. I'll break up the day by joining one of the sports clubs, going for a jog, or just grabbing dinner somewhere by the office. After the break, I'll get back to work, and I find that I'm quite productive. It's a great time to get in a four or five hour code sprint, where I can usually create a rough version of my newest idea.

Last week, I put together a new code for managing Fortran namelists (see an example of a namelist [here](http://jules-lsm.github.io/vn4.2/namelists/intro.html)), which is available as
[NameListManager](https://github.com/william-dawson/NameListManager) on github. Fortran lacks a standard library with  many built in features that users of other languages have come to expect, but surprisingly it comes with a very easy to use I/O routines for a specific type of file.

My group's code uses these files to specify the assortment of parameters needed by our quantum chemistry program. However, I noticed a real problem with this approach. Each developer would add new parameters to the code, and have to hardcode the default values inside the reader routines. This meant you had to look at the source code every time you wanted to know how an input parameter works.

Of course, documenting the code with doxygen would help. You would automatically get a website with a description of all the variables according to your in code comments. But consider this stub program for reading a namelist:

```Fortran
MODULE InModule
 IMPLICIT NONE
 REAL(8) :: converge_tol !< The maximum iterations to perform.
 INTEGER :: max_steps !< Tolerance to determine if a calculation converged.
CONTAINS
 SUBROUTINE ReadIn
   INTEGER, PARAMETER :: IO = 1
   NAMELIST /INPUT/ converge_tol, max_steps
   converge_tol = 1.0D-10
   max_steps = 100
   OPEN(UNIT=IO, FILE="INPUT", STATUS='OLD', ACCESS='SEQUENTIAL')
   READ(IO, INPUT)
   CLOSE(IO)
 END SUBROUTINE
END MODULE InModule
```

Now naturally we have documented the variables `converge_tol` and `max_steps` by writing some doxygen style comments beside them. But notice the problem with this approach. The initialization is separated from the documentation, which means you might end up with inconsistent default values. If you set the default values at the start:

```fortran
REAL(8) :: converge_tol = 1.0D-10
```

Subsequent calls to `ReadIn` won't reinitialize the default value. So naturally I overreacted to this problem and decided to build an entire package to manage this. `NameListManager` solves this problem in two steps. First, you begin by writing an XML file which describes the valid input parameters.

```xml
<?xml version="1.0" encoding="utf-8" ?>
<input name="example" lang="en" mpi="false">
  <group name="input">
    <description_list>
      <description lang="en">
        A module for storing the input parameters.
      </description>
    </description_list>
    <element_list>
      <element name="max_steps">
        <datatype>integer</datatype>
        <description_list>
          <description lang="en">
            The maximum iterations to perform
          </description>
        </description_list>
        <default>100</default>
      </element>
      <element name="converge_tol">
        <datatype>real(8)</datatype>
        <description_list>
          <description lang="en">
            Tolerance to determine if a calculation converged.
          </description>
        </description_list>
        <default>1.0D-10</default>
      </element>
    </element_list>
  </group>
</input>
```

Then if if you run the namelist manager on this input file, it automatically generates the source code for storing, initializing, and reading in these parameters. It also generates restructed text documentation that you might add to a sphinx project. Now when you want to add new input parameters to your program, you don't even have to look at the source code, but instead you spend time describing your new parameter.

One final note: why did I pick XML for input format? Of course, XML has the downside of being verbose, and isn't as trendy as say JSON or YAML. But (as described in earlier blog posts...) XML has the benefit of a schema file. I've included in `namelistmanager` a xml schema, which your input file is automatically checked against.

No doubt there are other great ways to accomplish this task, but if this light weight approach looks good to you, please clone this project on github. I hope to add this program to something [PyPI](https://pypi.org/) as well in the near future.
