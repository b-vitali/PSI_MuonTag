/// \file
/// \ingroup Tutorials
/// Example of `rootlogon.C`.
/// The macro `rootlogon.C` in the current working directory, is executed when
/// `root` starts unless the option `-n` is used.
///
/// \macro_code
///
/// \author Rene Brun

{
   printf("You again eh ? Keep it up...");
   printf("\nHopefully it will work this time around !\n\n");

   gStyle->SetPadGridX(kTRUE);
   gStyle->SetPadGridY(kTRUE);

}

