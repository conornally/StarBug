# Todo

work on todo list

main file:
> daophot wrap
> psf build
> prescience //or does this go into a file of its own 

alsclass.py:
> need to make sure that if i was to do any python based photometry, it would be able to handle not being taken from a file
> maybe a subclass?

fileio.py:
this is a os wrapper around the classes, ideally the classes wont have to know where to put things themselves //keep it general
> pass in filenames, pass back fits/als objects?
> more complex.. paths from files, paths from config, glob paths

> output to ./out/
> ths will need to know current working dir, and save dir etc

common_tasks.py
> dark frame flat field etc
> 

If i wanted to run pre science on a list of files
> config option to read from file     //.1
> create fitsclass for each, conduct all the stuff
> output to *separate* file >> this wont be too bad cause of self.name

