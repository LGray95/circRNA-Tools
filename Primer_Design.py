# Setting Up reverse complementary function

def reverse_complementary(inseq):
    outseq = inseq.replace("A","t").replace("T","a").replace("C", "g").replace("G", "c").upper()[::-1]
    return outseq

# Set up Primer Design function
def Primer_Design(circRNA_seq):

    # Count 10 nucleotides from start of sequence
    five_end = circRNA_seq[0:10]
    # Count 10 nucleotides fro the end of sequence
    three_end = circRNA_seq[-10:]


    # Creating Reverse Primer
    Reverse_Primer = reverse_complementary(three_end + five_end)

    # Finding the index 200 nucleotides from the end of circRNA
    index = (len(circRNA_seq)-180)
    Forward_Primer = circRNA_seq[index:(index+20)]

    # Printing output to user
    print("Forward Primer: 5'-",Forward_Primer,"-3'")
    print("Reverse Primer: 5'-",Reverse_Primer,"-3'")
    print("(Reverse Complementary to 5'-",three_end+five_end,"-3')")
    print("5'-", five_end)
    print("3'-", three_end)



Primer_Design("input_sequence")
