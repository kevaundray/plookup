## README 

POC of https://github.com/AztecProtocol/plonk-with-lookups/blob/master/PLONK-with-lookups.pdf

_Any misinformation presented here is the sole responsibility of me not understanding the protocol and does not reflect the protocol itself._

The base protocol is a multiset equality argument for showing that a multiset _f_ is contained in _t_.

_Use At Your Own Risk_

## Example 

```rust

// First instantiate a Lookup table. In our case, we have a 4 bit XOR Table
// However, the paper generalises to any multivariate function.
let table = XOR4BitTable::new();

// Now pass the table to the Lookup struct which will build the witness accordingly depending on your reads
let mut lookup = LookUp::new(table);

// Read will first check if (16 XOR 6) is in the 4-bit XOR table
// If it is, it will fetch the output and add the input and output to the witness
// A boolean is returned to indicate whether the input combination existed
// In the below example, since 16 cannot be represented using 4-bits, the witness would not have changed.
let added = lookup.read(&Fr::from(16), Fr::from(6)));

// Since 8 XOR 10 is available in the 4-bit XOR table
// 8, 10 and 8 XOR 10 will be added to the witness
lookup.read((Fr::from(8), Fr::from(10)));

// Alternatively, one can add the witness values directly without checking the table.
// Since there is a check that Z(X) was created correctly, this will fail on the prover side, if the values added are inconsistent with the table.
lookup.left_wires.push(Fr::from(1000));
lookup.right_wires.push(Fr::from(889));
lookup.output_wires.push(Fr::from(1234));

// Once all reads have been made. You can create a proof, that all of the witness values are indeed
// in the table. We create a random challenge by using a transcript object.
let mut transcript = Transcript::new(b"lookup");
let proof = lookup.prove(&proving_key, &mut transcript); 

```

## How does this integrate into PLONK?

- Let's use M(X) to denote the multiset equality polynomial. In the above paper, this is represented using Z(X), however in PLONK Z(X) is the accumulator for the permutation polynomial.
- Lets also uses Q(X) for the quotient polynomial as in PLONK because t(X) has a different meaning in the above paper.


- 3 commitments are added to the proof size: h_1(X), h_2(X), M(X)
- 2 Field elements are added to the proof size: evaluations of \hat{M(X)} and \hat{M(Xg)}
- The prover would need to modify the linearisation polynomial to account for the extra terms in the quotient polynomial
- The prover would need to modify the opening polynomial to include the openings for M(X) and M(Xg)
- The quotient polynomial Q(x) is extended to check the necessary items in this protocol. It is orthogonal to the custom gates.


## Caveats

- The Quotient polynomial is not split into degree-n polynomials, so the SRS is not linear in the number of reads. This can be fixed quite easily in the POC, by doing the same technique that PLONK did.

- In the protocol, we use a random challenge `alpha` to _fold_ the lookup table into a vector. In an interactive setting, this is fine as `alpha` is random. In a Non-interative setting, this becomes a problem, as the transcript is empty, prior to generating `alpha`. The consequence of this is that one will need to embed this protocol in another protocol, so that you have sufficient entropy to generate alpha.

- This POC does not aggregate the witness. We will therefore have a witness per commitment, in reality, we will have one witness for all polynomials, since they are all evaluated at the same point. This was done for easier debugging. 