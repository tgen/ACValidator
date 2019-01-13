Instructions for testing ACValidator using the test data provided
-----------------------------------------------------------------

Test data are available here: http://tools.tgen.org/Files/ACValidator_test_data/ 
An input SAM file from one of our simulations is provided here for testing: http://tools.tgen.org/Files/ACValidator_test_data/InputFile/ 
Coordinates that can be used for testing are in http://tools.tgen.org/Files/ACValidator_test_data/test_coordinates.txt

1) To run directly from source code
a) For a single coordinate:
python ACValidator.py -i positive1 -c 8:103312227-103372418 -w 300 --log-filename Log.txt

b) For all the coordinates in test_coordinates.txt, using the launcher: 
(make sure ACV_launcher.sh and ACValidator.py are in the same folder if running directly from source code)
./ACV_launcher.sh positive1 test_coordinates.txt window_size Log.txt


2) To run after pip install
a) For a single coordinate:
ACValidator -i positive1 -c 8:103312227-103372418 -w 300 --log-filename Log.txt

b) For all the coordinates in test_coordinates.txt:
for coordinate in `cat test_coordinates.txt`; do ACValidator -i positive1 -c ${coordinate} -w 300 --log-filename Log.txt; done 
