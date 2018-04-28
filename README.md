# Asset-Backed-Security-ABS-
Pricing Asset Backed Securities (ABS)

Consider a 30-year MBS with a fixed 𝑊𝐴𝐶 =8% (monthly cash flows starting in January of this year). The Notional Amount of the Loan is $100,000. Use the CIR model of interest rates 𝑑𝑟𝑡=𝜅(𝑟̅−𝑟𝑡)𝑑𝑡+𝜎√𝑟𝑡𝑑𝑊𝑡 with 𝑟0=0.078,𝑘=0.6, 𝑟̅=0.08,𝜎 =0.12

1. Consider the Numerix Prepayment Model.

(a) Compute the price of the MBS using this model for prepayments. The code should be generic: the user is prompted for inputs and the program runs and gives the output.

(b) Compute the price of the MBS for the following ranges of the parameters: 𝑘 in 0.3 to 0.9 (in increments of 0.1) and draw the graph of the price vs. 𝑘.

(c) Compute the price of the MBS for the following ranges of the parameters: 𝑟̅ in 0.03 to 0.09 (in increments of 0.01) and draw the graph of the price vs. 𝑟̅.

2. Consider the PSA Model of prepayments.

(a) Compute the price of the MBS using the PSA model for Prepayments. The code should be generic: the user is prompted for inputs and the program runs and gives the output.

(b) Compute the price of the MBS for the following ranges of the parameters: 𝑘 in 0.3 to 0.9 (in increments of 0.1) and draw the graph of the price vs. 𝑘.

3. Compute the Option-Adjusted-Spread (OAS) for the Numerix-Prepayment model case with the Market Price of MBS being $110,000.

4. Compute the OAS-adjusted Duration and Convexity of the MBS, considered in the previous question.

5. Consider the MBS described above and the IO and PO tranches. Use the Numerix-Prepayment Model and price the IO and PO tranches for: 𝑟̅ in 0.03 to 0.09 range, in increments of 0.01.
