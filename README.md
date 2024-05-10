# Finch-to-Finch Symmetric Compiler 
Taking advantage of symmetry in matrices saves a factor of two, but taking advantage of symmetry in a tensor of order $n$ can save a factor of $n!$ in memory accesses and operations. With this compiler, you can capitalize on symmetry in tensor kernels, resulting in savings ranging from 1.3x for SpMV to 7.8x for a 4-dimensional MTTKRP kernel!

## Functionality 
This compiler utilizes RewriteTools, a utility for term rewriting that provides a language for finding subexpressions that satisfy specific conditions and applying predefined
transformations on the matches. The compiler operates in two phases:

1. _Symmetrization_: a symmetric kernel that only accesses the canonical triangle of the symmetric inputs is generated (saving memory bandwidth). 
2. _Optimization_: a series of transforms are performed on the symmetric kernel to filter redundant computations (saving computational bandwidth).

The compiler takes in a Finch program expression as input and returns an executable Finch program as output. 

## Usage
You can load and call the compiler from the Julia command-line as follows: 
~~~
julia> include("PATH/TO/REPO/compiler/symmetrize.jl")
~~~

The entry-point to the compiler is a function `symmetrize` that takes in three inputs: an assignment operation written as a Finch program, a list of the symmetric tensor inputs, and a list of indices ordered in the desired loop order (outermost to innermost). 

The following is an example for the SpMV kernel: 
~~~
julia> y = :y; A = :A; x = :x;

julia> i = index(:i); j = index(:j);

julia> ex = @finch_program y[i] += A[i, j] * x[j];

julia> symmetrize(ex, [A], [i, j])
Finch program: for j = virtual(Dimensionless), i = virtual(Dimensionless)
  if <=(i, j)
    let x_j = x[j], A_ij = A[i, j]
      begin
        if !=(i, j)
          begin
            y[i] <<+>>= *(A_ij, x_j)
            y[j] <<+>>= *(A_ij, x[i])
          end
        end
        if ==(i, j)
          begin
            y[i] <<+>>= *(A_ij, x_j)
          end
        end
      end
    end
  end
end
~~~


