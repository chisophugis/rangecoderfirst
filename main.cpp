// Range decoding
//
//
// The basic idea of range coding:
// ===============================
//
// The basic idea of range coding is to represent a message as a binary
// fraction in the range [0,1). Suppose you are given a discrete CDF with
// stops [0, 3/4, 1), and a message corresponding to the value V = 0.1111
// (all such fractions will be binary in this document).
// To decode this message, you look at V, and compare it with the CDF. If V
// lies in the range [0,3/4=0.11), then you decode a symbol 'a', if it lies
// in the range [3/4=0.11,1), then you decode 'b'.
// In this case, 0.11 <= V = 0.1111 < 1, so you decode 'b'.
// There are two interpretations of what you do next:
//   1. You "scale" the range [3/4,1) back to the range [0,1), and apply
//      the same scaling to V. This has a nice recursive flavor since now
//      you are again set up with [0,1) and another value V.
//      - This interpretation is generally useful from a low-level
//        implementation perspective, since you can always keep the
//        "tangible" range [0,1), which is roughly what you represent in
//        your data structure.
//   2. You "zoom in" on [3/4,1), with V lying within that range.
//      - This interpretation is generally useful from a theoretical
//        perspective, but it adapts itself to describing the artifacts of
//        a fixed-precision implementation better, since "zooming in" on a
//        long stream of bits is basically equivalent to left
//        shifts/incrementing pointers, which doesn't lose precision.
// Both interpretations are equivalent, and both are useful ways of
// thinking about range coding. Ultimately the low-level code for the
// implementation looks sort of like a demented "scale" operation.
//
// From the "scale" perspective, your next step to decode this message is
// to scale the range [3/4,1) back to [0,1), which involves subtracting
// 3/4=0.11 and scaling back up by the reciprocal of 1-3/4=1/4, i.e.
// multiplying by 4.
// This operation does the following to V:
//   0.0011 = V - 0.11
//   0.11 = 0.0011 * 4.
// So the new value of V is 0.11. You can probably already guess that the
// next decoded symbol is 'b' again (if not, notice that 3/4 <= 0.11 < 1).
// Typically, there is some contextual way for the decoder to know when to
// stop; the decoder effectively goes through the motions of a very tight
// synchronized ballet dance with the faraway encoder (see the next
// paragraph), so this is not very difficult to arrange.
//
// Many of the steps above required
// 1. arbitrary precision arithmetic, and
// 2. multiplications/divisions.
// We gain speed by using alternative operations that operate on fixed
// precision numbers and which avoid multiplications/division (except where
// they correspond to bit shifts).
// Additionally, in a "real" range coder the CDF that is used can change
// from symbol to symbol; the only caveat is that both the encoder and
// decoder must use the exact same CDF's in the exact same order.
// In general, the decoder must do everything in *exactly* the same way as
// the encoder, including the *exact* same "errors" due to the
// fixed-precision arithmetic; everything must be identical, or else things
// will go catastrophically wrong.
//
// The general model:
// ==================
//
// Notionally, we are given a long (but finite) string of bits, which is
// interpreted as a base-2 fraction V = 0.vvv...vvv00... (where the v's
// represent arbitrary bit values). The name `V` just stands for "value".
//
// Obviously we do not want to do arithmetic on arbitrarily long strings of
// bits, so instead we operate on a finite-precision "subview" of V (I
// would like to use the word "window", but that already has a slightly
// different meaning in the Daala range coder).
// In the case of this range coder, modeled after the Daala range coder,
// this subview is 16 bits wide (but the description given here generalizes
// easily), and is represented as the following:
//
// R = 1rrr rrrr rrrr rrrr (i.e. `rng` in the Daala source)
// A = aaaa aaaa aaaa aaaa (i.e. `dif >> 16` in the Daala source)
//   satisfying the following invariants:
//     * 0 <= A < R
//     * 2^15 <= R < 2^16
//     * The high bit of R is 1 (i.e. R >= 2^15). This condition guarantees
//       that A can take on at least 2^(16-1) values, which provides a
//       precision guarantee.
//     * Whenever you shift A to the left, you must shift in corresponding
//       bits from V (shifting A to the left is done during
//       "renormalization").
//       In the terminology developed below, this consdition is equivalent
//       to: V - Z = InfPrec(A) + epsilon, where epsilon < InfPrec(1).
// Graphically:
//   [--------|---------------)
//   0        A        2^15 < R < 2^16
//
// This range [0,R) corresponds to an interval on the (abitrary-precision)
// real interval [0,1). Throughout the decoding process, this
// correspondence changes, as will be exhibited in the diagrams throughout.
// We will use the notation InfPrec(X), where X is a fixed-precision number
// like R or A, to represent the corresponding number in the
// arbitrary-precision interval [0,1).
//
// The initial state of the decoder is such that InfPrec(R) = 1.0, that is,
// the finite-precision interval [0,R) corresponds to the
// arbitrary-precision real interval [0.00...,1.00...):
//         R  = 1 000 0000 0000 0000
// InfPrec(R) = 1.000 0000 0000 0000 = R * 2^(-(16-1))
//         A  = 0 vvv vvvv vvvv vvvv
// InfPrec(A) = 0.vvv vvvv vvvv vvvv = A * 2^(-(16-1))
//         V  = 0.vvv vvvv vvvv vvvvvvvvv...
//         Z  = 0.000 0000 0000 000000000... = InfPrec(0)
//              ^
//              `-- k = 0 and points here.
//              NOTE: The grouping of digits is purely cosmetic, and does not
//              imply any sort of "alignment". e.g. in `0.000...00 xxxx`, the
//              number of leading 0's is arbitrary, and not required to be a
//              multiple of 4.
// The number Z satisfies the relation "V - Z = InfPrec(A) + epsilon",
// where "epsilon < InfPrec(1)", which basically means "the only nonzero
// bits of epsilon are to the right of our 16-bit fixed-precision window",
// as represented in this diagram of the "general" state of the decoder
// when all the invariants are satisfied:
//         R  =            1rrr rrrr rrrr rrrr
// InfPrec(R) = 0.000...00 1rrr rrrr rrrr rrrr
//         A  =            aaaa aaaa aaaa aaaa
// InfPrec(A) = 0.000...00 aaaa aaaa aaaa aaaa 000...
//         V  = 0.vvv...vv vvvv vvvv vvvv vvvv vvv...
//         Z  = 0.vvv...vv zzzz zzzz zzzz zzzz 000... = InfPrec(0)
// InfPrec(1) = 0.000...00 0000 0000 0000 0001 000...
//   epsilon  = 0.000...00 0000 0000 0000 0000 vvv...
// In general, the bits z in the above diagram are nonzero (see the
// description below of how A and R are updated to understand why; consider
// s_K == 1 and A == aaaa...0001).
// Note that InfPrec([0,R)) = [Z, Z + InfPrec(R)) represents a subrange of
// the infinite precision range [0,1).
//
// Decoding a range:
// -----------------
//
// Now, the fundamental operation on this subview consists of "decoding a
// range", as follows:
//
//   Input: A set of distinct stops s_i satisfying:
//     s_i < s_{i+1} and s_0 = 0, s_N = R
//   These are to be interpreted as the CDF of a probability distribution.
//   Graphically, the input is like this:
//     [----------|--------|--------...---|------------)
//     0=s_0      s_1      s_2      ...   s_{N-1}      R=s_N
//   (note that all of the intervals [s_i,s_{i+1}) are disjoint and
//   nonempty).
//   Then there is a unique K such that [s_K,s_{K+1}) contains A.
//   >>> The decoded value is K, which satisfies 0 <= K < N. <<<
//   Now, as far as the user of the decoder is concerned, this is all that
//   matters, since their symbol K has been decoded. However, the decoder
//   must update its state, with a process called "renormalization".
//
// Renormalization:
// ----------------
//
//   The first step of updating R and A is as follows:
//     R' = s_{K+1} - s_K
//     A' = A - s_K
//     Note that we still have A' < R' since
//         A' < R' <==> A' + s_K < R' + s_K <==> s_{K+1} < A"
//     Graphically:
//     [--------...--------|------------|--------------...---------)
//     0=s_0    ...        s_K          s_{K+1}        ...         R
//                         [------------)
//                         0'     A'    R' = s_{K+1} - s_K
//   However, this update may break the invariant set forth above that the
//   high bit of R must be 1. A typical case may be:
//         R  =            1rrr rrrr rrrr rrrr
//         R' =            0001 ssss ssss ssss (Unrelated to the s_i. Oops!)
//         A  =            aaaa aaaa aaaa aaaa
//         A' =            000b bbbb bbbb bbbb
//         V  = 0.vvv...vv vvvv vvvv vvvv vvvv vvv...
//         Z  = 0.vvv...vv zzzz zzzz zzzz zzzz 000... = InfPrec(0)
//         Z' = 0.vvv...vv vvvz zzzz zzzz zzzz 000... = Z + InfPrec(s_K)
// InfPrec(1) = 0.000...00 0000 0000 0000 0001 000...
//   epsilon  = 0.000...00 0000 0000 0000 0000 vvv...
//                         ^
//                         `-- k points here.
//                         Note that the leading bit of R' is not 1! The
//                         invariants do not hold!
//  To restore the invariants, we must:
//    * shift `R'` to the left to form R~,              Side note:       ~  ~  ~
//    * "shift in bits" from V into A' to form A~       typographically: R, A, k
//    * add the shift amount to k to form k~
//  (All these X~ will be come just X for the next decoding step).
//  The results is that (say there are 3 leading zeros in R'):
//         R  =            1rrr rrrr rrrr rrrr
//         R' =            0001 ssss ssss ssss
//         R~ =               1 ssss ssss ssss 000
//         A  =            aaaa aaaa aaaa aaaa
//         A' =            000b bbbb bbbb bbbb
//         A~ =               b bbbb bbbb bbbb vvv
//         V  = 0.vvv...vv vvvv vvvv vvvv vvvv vvv...
//         Z  = 0.vvv...vv zzzz zzzz zzzz zzzz 000... = InfPrec(0)
//         Z' = 0.vvv...vv vvvz zzzz zzzz zzzz 000... = Z + InfPrec(s_K)
//         Z~ = 0.vvv...vv vvvz zzzz zzzz zzzz 000... = Z' (Why?)
// InfPrec(1) = 0.000...00 0000 0000 0000 0001 000...
//   epsilon  = 0.000...00 0000 0000 0000 0000 vvvvvv... < InfPrec(1)
// InfPrec(1~)= 0.000...00 0000 0000 0000 0000 001000...
//   epsilon~ = 0.000...00 0000 0000 0000 0000 000vvv... < InfPrec(1~)
//                         ^  ^
//                         |  `-- k~ points here.
//                         `------ k points here.
//            NOTE: Z' = Z~ because the low bits of A correspond exactly
//            with the corresponding bits of V.
//   The "tilde" values now satisfy our invariants and we set:
//     R := R~
//     A := A~
//     There is no variable `k` in the program, but notionally:
//     k := k~ (this changes the meaning of InfPrec!)
//   Basically, we make sure that the condition "the highest bit of R is 1"
//   is true by redefining what InfPrec means, which basically corresponds
//   to shifting our subview further to the right (or shifting the values
//   of R~ and A~ to the left).
//
//  TODO: to avoid "destructive assignments" in the math, could use have a
//  subscript "n" on everything indicating "that variable, after n symbols
//  have been decoded". That might clutter things up. Definitely not worth
//  it in the plaintext version, but might be fine in TeX. That would avoid
//  the "tilde" variables.
//
// Reinterpreting this description into a practical implementation:
// ================================================================
//
// The above description focused a lot on the relationship of the
// fixed-precision values R and A with notional infinite precision
// counterparts. In this section, we will not use InfPrec even once!
//
// The input to the range coder is a buffer of bytes, called V (this is a
// monospace "code" V, not the "math" variable V in the preceding section).
// V is as follows in memory:
//
// |vvvv vvvv|vvvv vvvv|vvvv vvvv|vvvv vvvv|...
// ^
// |
// &V
// The issue of bit-endianness comes up. If you zoom into each of the bytes
// above, it is numbered as follows:
//    |0 1 2 3  4 5 6 7|
// This turns out to be the "natural" ordering in the sense that it does
// not require bit reversals as bits are read from memory. E.g.
// * `(V[0] & 0b1000_0000) >> 7` is the first bit in the bitstream,
// * `(V[0] & 0b0100_0000) >> 6` the second,
// * etc.
// (Note: reading the bits like that would be quite naive)
//
// We use a type `Integer` to represent the types of the values R and A.
// The type `Symbol` represents an integer holding the value of a decoded
// symbol.
//
// The decoding procedure is as follows:
//
// * Subview represents a class with two fields, R and A, and whose
//   instances *always* satisfy the invariants on R and A. It is passed
//   around by value. If R and A are each 16 bits, then Subview will occupy
//   32 bits.
// * Lo and Hi correspond to s_K and s_{K+1} in the description of the
//   previous section. The process of computing these values in `findSymbol`
//   (i.e. searching the CDF) needs to be discussed later, since in general
//   the CDF is not provided as the values s_i representing fractional
//   stops in the range [0,R), and so a "range expansion" process is needed
//   to scale them to this range (naively, the range expansion could be
//   just an appropriate multiplication/division, but for performance
//   reasons a different expansion technique is used).
// * Note that the use of a CDF below can be optimized in some cases. This
//   is the "most general" representation of the algorithm, but there are
//   important special cases that can be optimized.
// * BitstreamReader corresponds to a class that reads a bitstream bits at
//   a time; this operation can be optimized and will be discussed later.
//   Only BitstreamReader reads from V.
//
// pair<Symbol, SubView> decodeRange(SubView SV, BitstreamReader &BR, CDF &C) {
//   LookupResult Result = findSymbol(CDF, SV);
//   return {Result.Sym, renormalizeDec(SV, BR, Result.Lo, Result.Hi)};
// }
//
// The renormalization procedure follows immediately from the definitions
// of R', A', R~, and A~:
//
// SubView renormalizeDec(SubView SV, BitstreamReader &BR,
//                        Integer Lo, Integer Hi) {
//   Integer RPrime = Hi - Lo;
//   Integer APrime = SV.A - Lo;
//   Integer LZCNT = lzcntNonZero(RPrime);
//   // Note that BR.readBits(0) should be a nop.
//   Integer RTilde = (RPrime << LZCNT);
//   Integer ATilde = (APrime << LZCNT) | BR.readNBits(LZCNT);
//   return SubView{RTilde, ATilde};
// }
//
// (Implementation note: the core of the preceding function is essentially
// 5 x86 instructions "SUB; SUB; BSR; SHL; SHLD;")
//
//
// "Range expansion", i.e. scaling symbol probabilities to the range [0,R):
// ------------------------------------------------------------------------
//
// The values Lo and Hi in the description above are the numerators of
// fractions with denominator R, representing a symbol probability.
// This is necessary for the decoding procedure to operate correctly.
// However, note that R changes throughout the decoding procedure.
// This conflicts with the fact that symbol probabilities are typically
// expressed as fixed point fractions with a fixed denominator (e.g. as Q15
// fixed point numbers), or a denominator whose value is independent of R
// (e.g. the numerators may be frequency counts, and the denominator the
// total frequency count).
//
// This issue is addressed with a "range expansion" operation, which
// expands a fraction P/Q representing a symbol probability with a fixed
// denominator into a fraction P'/R that in some sense "corresponds" to
// P/Q. Naively, we would choose P' so that
//  P     P'
// --- = ---
//  Q     R
// and do the appropriate multiplications and division to derive a value of
// P' from P (i.e., multiply by "R/Q") that is fairly close to the ideal
// value (if P/Q is incommensurable with 1/R, then it may not be exact).
// However, this is fairly expensive; we would like to get away with doing
// no multiplications or divisions (unless they can be implemented as bit
// shifts). The only paper where I have seen the issue of efficient range
// expansion addressed is [SM98] and
// <http://tools.ietf.org/html/draft-terriberry-codingtools-00> which
// references [SM98].
//
// [SM98]     Stuiver, L. and A. Moffat, "Piecewise Integer Mapping for
//            Arithmetic Coding", Proc. of the 17th IEEE Data
//            Compression Conference (DCC'98) pp. 1--10, March/
//            April 1998.
//
// Daala uses a "piecewise linear" mapping similar to the one described in
// [SM98], but "reversed" in the sense that a prefix of the range is
// expanded by 2x and a suffix is expanded 1x ([SM98] does 1x on a
// prefix and 2x on a suffix).
// The denominator Q is assumed to satisfy Q <= R <= 2Q; this can be
// accomplished by shifting Q to the left appropriately (and shifting the
// corresponding P values as well), as long as Q <= 2^15.
// In fact, a value Q = 2^15 always satisfies this since one of the
// invariants on R is 2^15 <= R < 2^16.
// The range expansion operation is given by:
//   P |--> P + min(x, R - Q)
// In other words, P/Q is turned into
//  P        P + min(P, R - Q)
// --- |--> -------------------
//  Q              R
//
//  This corresponds to a piecewise linear mapping that can be graphically
//  interpreted as follows:
//    d = 2Q - R
//    R - Q = Q - (2Q - R) = Q - d
//
//    0       Q-d    Q            2Q
//    [--------)-----)-------------)
//    |         \     \            .
//    |          \     \           .
//    |  2x       \     \          .
//    |  |    ^    \ 1x  \         .
//    |  V    |     \     \        .
//    |     (1/2)x   \     \       .
//    |               \     \      .
//    [----------------)-----)     .
//    0                R-d   R     .
//                     [-----)-----)
//                     0     d     2*d
// To get a feel for this mapping, observe that when Q=2^15 and R=2^15,
// then d = 2Q - R = 2^15, and the entire mapping is a "1x" expansion,
// i.e. and identity mapping.
// On the other hand, when R=2^15 and Q=2^14 (the minimum allowable Q for
// such an R), then d = 2Q - R = 0 and the entire mapping is a "2x"
// expansion.
//
// It is often useful to have the inverse operation so that computations
// can be done in the "P/Q" space.
// For example, A/R can be mapped into a fraction P/Q with one computation,
// and this is used to search a CDF with denominator Q rather than
// converting each value in the CDF to the form P'/R.
// Such an inverse mapping is:
//    x |-> max(x/2, x-(R-Q)) = max(x/2, x-(Q-d))
// Where d is as in the figure above (this equation uses the variable `d`
// from the forward mapping, which diverges slightly from the variable
// named `d` in e.g. `od_ec_decode_cdf`, which represents R-Q)
// This mapping consists of region of slope 1/2 until a "crossover point"
// when x/2 = x-(Q-d), which resolves to x = 2(Q-d) = R-d after which the
// mapping has slope 1. This corresponds to going from the bottom line to
// the top line in the above figure of the forward mapping, where the 2x
// expansion is now a (1/2)x contraction.
//
// Reading from the bitstream:
// ---------------------------
//
// Reading a bitstream a bit at a time is a simple exercise that probably
// every programmer has done at some point.
// The approach described here is probably the second thing that any
// programmer would come up with (after the really simple approach of
// maintaining a "byte index" and a "bit index").
//
// The idea is to have a single word Buf, of a size natural for the
// implementation architecture (32 or 64 bits, most likely), that acts as a
// "bit buffer". A memory read from the actual input byte stream V is done
// in word-aligned reads that refill Buf entirely.
//
// Recall from above that the interface to a BitstreamReader is
// `Integer readNBits(int N)`, where N is the number of leading zeros of
// R' = Hi - Lo.
// Note that a range [Lo,Hi) with Lo == Hi has probability
// zero and cannot be decoded, so it follows that R' >= 1 [*].
// It follows that for a 16-bit subview, N will be at most 15 bits.
// If Buf contains the requested number of bits, then they are shifted out
// and returned (updating Buf in the process of shifting the bits out).
// Otherwise, any bits available are shifted into a temporary (thus leaving
// Buf = 0), and a single aligned read of sizeof(Buf) bytes completely
// refills buf from the byte stream V, from which the remaining deficit of
// bits (i.e., N - # bits available in the temporary) is fulfilled in the
// usual way.
// (There is also the infrequent (once per packet) "end of stream"
// condition, where a read of the word size would read past the end of the
// array V, which must be guarded against, but an implementation should
// avoid inducing a large cost for this range check to the detriment of
// decoding performance in the general case. The "raw bits" encoded at the
// end of a Daala packet (a topic for another section) might naturally
// provide sufficient padding to guarantee this in all cases of practical
// interest).
//
// This can also be interpreted as a FIFO of bits, where the first 16 bits
// of the FIFO represent A, and forming A' from A corresponds to
// subtracting Lo from the first 16 bits of the FIFO, and then "shifting A'
// to the left" to form A~ corresponds to dropping the leading 0 bits from
// the FIFO.
// This interpretation is likely to be useful for hardware implementations.
// In this scheme R will have to be transported alongside the FIFO.
// The notional number "epsilon" in the theoretical description above
// corresponds to the "tail" of this FIFO (all but the first 16 bits).
//
// (Note that Daala's current implementation effectively keeps a single
// 32-bit word "dif = A:Buf" which corresponds to a Buf of only 16 bits.
// Daala's current implementation also refills the buffer one byte at a
// time rather than using larger aligned loads).
//
// [*] In the theoretical discussion above, we required that the ranges
// [s_i,s_{i+1}) be nonempty. In practice, these ranges just can be left
// empty, corresponding to a probability of zero for the corresponding
// symbol, if that is convenient (such as disabling a feature without
// having to change the decoding code).
//
// General Note: everything here depends *crucially* on the encoder doing the
// range expansion and renormalization in *exactly* the same way.
//

// Feedback:
// Biggest point of confusion is InfPrec (which I had defined wrong! it is
// not just a multiplication by 2^(-(k-1))).
// Also Z, which corresponds to what in the literature is called L.
