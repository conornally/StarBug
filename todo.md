# Todo

work on todo list

main file:
> daophot wrap
> psf build
> more group interaction
>> done d, move

fitsclass.py
> shift file i/o

alsclass.py:
> need to make sure that if i was to do any python based photometry, it would be able to handle not being taken from a file
> maybe a subclass?
> als2reg

fileio.py:
this is a os wrapper around the classes, ideally the classes wont have to know where to put things themselves //keep it general

common_tasks.py

psf selector 
> (works for both daophot and photutils?) im sure i can get that to work ahaha...






#################################################
If i wanted to run pre science on a list of files
> config option to read from file     //.1
> create fitsclass for each, conduct all the stuff
> output to *separate* file >> this wont be too bad cause of self.name

