## README 

POC of https://github.com/AztecProtocol/plonk-with-lookups/blob/master/PLONK-with-lookups.pdf

_Any misinformation presented here is the sole responsibility of me not understanding the protocol and does not reflect the protocol itself._

The base protocol is a multiset equality argument for showing that a multiset _f_ is contained in _t_.


## How does this integrate into PLONK?

- Let's use M(X) to denote the multiset equality polynomial. In the above paper, this is represented using Z(X), however in PLONK Z(X) is the accumulator for the permutation polynomial.
- Lets also uses Q(X) for the quotient polynomial as in PLONK because t(X) has a different meaning in the above paper.


- 3 commitments are added to the proof size: h_1(X), h_2(X), M(X)
- 2 Field elements are added to the proof size: evaluations of \hat{M(X)} and \hat{M(Xg)}
- The prover would need to modify the linearisation polynomial to account for the extra terms in the quotient polynomial
- The prover would need to modify the opening polynomial to include the openings for M(X) and M(Xg)
- The quotient polynomial Q(x) is extended to check the necessary items in this protocol. It is orthogonal to the custom gates.