# Todo

>catclass load sb files // done basic
>need to fix sourceclass calling nan error = 0


work on todo list

main file:
> daophot wrap
> psf build
> more group interaction
>> done d, move
> tab complete in get_group and load

fitsclass.py
> shift file i/o

alsclass.py:
> need to make sure that if i was to do any python based photometry, it would be able to handle not being taken from a file
> maybe a subclass?
> scrap (for now all catalog exporting compression) just export the whole thing

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

