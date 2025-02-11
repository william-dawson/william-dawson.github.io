Title: Data Quality in the Fitting of Approximate Models: A Computational Chemistry Perspective
Date: 2025-02-11 11:00:00 +0900
Category: Lessons
Tags: publication
Summary: I share a paper we wrote about data quality in computational chemistry calculations, and spend some time hacking Libxc.

I had planned to write two blog post entries over the Christmas holiday, but ended up down for the count with the flu. Today is a holiday and it's snowing outside, so it seems like a chance to stay cozy inside and knock a post off the to do list. This entry is inspired by our recent publications entitled [Data Quality in the Fitting of Approximate Models: A Computational Chemistry Perspective](https://doi.org/10.1021/acs.jctc.4c01063). Our paper asks how important is it to have high quality data when fitting a new computational chemistry method. Surprisingly, the answer is: not as important as you'd think. This comes down to the balance between model limitations and overfitting; I'll leave the details to the paper itself.

Instead, for the blog entry I want to talk about hacking [Libxc](https://libxc.gitlab.io/). Libxc is one of the best libraries in our field, implementing the many density functional approximations out there in the literature. [It can be tricky to implement these approximations just by reading the papers
](https://doi.org/10.1063/5.0167763), so having a library to do all the work can save developers a lot of time. 

One of the functionals they implement is called [B97](https://doi.org/10.1063/1.475007), proposed by Axel D. Becke. B97 is a great example of an empirical density functional: thanks to Becke's linearization trick, you can expand the functional as a linear combination of multiple terms. This makes it very easy to fit to any dataset you have. Unsurprisingly, there are now many variants of the original B97. In Libxc, [there are currently 23 versions of it](https://github.com/ElectronicStructureLibrary/libxc/blob/master/src/gga_xc_b97.c). The source code here is very readable, so you can see the parameters for any given definition right away:
```
static const double b97_values[B97_N_PAR] =
  {0.8094, 0.5073, 0.7481, 0.0, 0.0,
   0.1737, 2.3487, -2.4868, 0.0, 0.0,
   0.9454, 0.7471, -4.5961, 0.0, 0.0,
   0.1943};
```
This values are identical to Table III in Becke's paper, except that we have the possibility for some extra values (there are 15 parameter and 9 parameter versions, with both possibily being hybrid functionals).

My only complaint about Libxc here is that these are const values! What if I want to modify them myself and create my own functional? This is where we start hacking. First, remove that pesky `const` keyword. Then add a function to allow us to modify the parameters:
```
void in_set_b97_values(const double* values, size_t size) {
    assert(size == B97_N_PAR); // Ensure size match
    for (size_t i = 0; i < size; ++i) {
        b97_values[i] = values[i];
    } 
} 
```
Now let's make this work with [PySCF](https://pyscf.org/). We need to get PySCF to do the proper wrapping of our C function into python, which can be done by modifying `pyscf/pyscf/lib/dft/libxc_itrf.c`:
```
void LIBXC_set_b97_values(const double* values, size_t size) {
        in_set_b97_values(values, size);
}
```
and then adding a function to `pyscf/pyscf/dft/libxc.py`:
```
def set_b97_values(v):
    n = len(v)    
    _itrf.LIBXC_set_b97_values((ctypes.c_double*n)(*v), ctypes.c_int(n))
```

That's all it takes to generate your own functional. This can be tested by verifying that we can reproduce another parameterization by modifying the values to match.
```
from pyscf import dft, gto
from pyscf.dft import libxc
from sys import argv

mol = gto.M(atom = 'O 0 0 0 ; '
                   'H 0.759 0 0.585 ; '
                   'H -0.759 0 0.585',
            basis = '6-31G*')
mol.verbose = 0
rks = dft.RKS(mol) 

if argv[1] == "b97":
    rks.xc = 'b97'
elif argv[1] == "b97_1":
    rks.xc = 'b97_1'
elif argv[1] == "fake":
    rks.xc = 'b97'
    libxc.set_b97_values([
       0.789518, 0.573805, 0.660975, 0.0, 0.0,
       0.0820011, 2.71681, -2.87103, 0.0, 0.0,
       0.955689, 0.788552, -5.47869, 0.0, 0.0,
       0.21
    ])
rks.kernel()
print(argv[1], rks.e_tot)
```
```
b97 -76.3802820297083
b97_1 -76.38208191613612
fake -76.38208191613609
```
These days it is probably difficult to simply reparameterize a functional and get a publication, for reasons detailed in our paper. But don't let that stop you from realizing your dream or having your own functional, especially when the code is this accessible. 

