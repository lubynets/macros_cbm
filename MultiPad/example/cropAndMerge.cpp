#include "MultiPicture.h"

// Compile this code using g++ cropAndMerge.cpp ~/path-to-MultiPad/MultiPicture.cpp -I ~/path-to-MultiPad
// And then run the executable ./a.out

int main() {

  MultiPicture mpic(3, 2);
  mpic.SetVerbose();

  // Set pictures in a grid
  // The first coordinate - column, the second - row
  // Start counting from 0 from the top left corner
  mpic.SetPictureName(0, 0, "topLeft.pdf");
  mpic.SetPictureName(0, 1, "bottomLeft.pdf");
  mpic.SetPictureName(1, 0, "topMid.pdf");
  mpic.SetPictureName(1, 1, "bottomMid.pdf");
  mpic.SetPictureName(2, 0, "topRight.pdf");
  mpic.SetPictureName(2, 1, "bottomRight.pdf");
  
  mpic.SetLeftMargins({0.14, 0.14, 0.14}); // Set margins as they are in root TCanvas
  mpic.SetTopMargins({0.07, 0.07}); // Set margins as they are in root TCanvas
  mpic.SetRightMargins(0.03); // if all margins are equal call function with a single argument
  mpic.SetBottomMargins(0.12); // if all margins are equal call function with a single argument

  mpic.SetLeftMarginShifts(0.03);
  mpic.SetBottomMarginShifts(0.01);

  // When using this functions, in the macro which produces input pdf or png
  // there must be explicit call of TCanvas::SetCanvasSize(1200, 800);
  mpic.SetPixelSizesX(1200);
  mpic.SetPixelSizesY(800);

  mpic.Run();

  return 0;
}
