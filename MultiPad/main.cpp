#include <iostream>
#include "MultiPicture.h"

int main() {

  MultiPicture mpic(3, 2);
  mpic.SetVerbose();

  mpic.SetPictureName(0, 0, "image1.png");
  mpic.SetPictureName(0, 1, "image2.png");
  mpic.SetPictureName(1, 0, "image3.png");
  mpic.SetPictureName(1, 1, "image4.png");
  mpic.SetPictureName(2, 0, "image5.png");
  mpic.SetPictureName(2, 1, "image6.png");

  mpic.SetLeftMargins(0.15);
  mpic.SetRightMargins(0.01);
  mpic.SetTopMargins(0.01);
  mpic.SetBottomMargins(0.08);

  mpic.Run();

  return 0;
}
