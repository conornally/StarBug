# StarBug

## Basic I/O

`> load`

Prompted with the file input command, this takes a variety of input styles. A single fits file can be loaded with the path to the file. Any number of files can be loaded with a list of paths, space separated. [Allow `*`.fits inputs]. A file containing list of paths can be input by specifying path to that file.
The next *input* is group name; this represents the name for the list (or single) fits object inside the program, and will be used later on for any analysis regime.

*NB:* There is special group names pre-loaded:

| Name | Used for              |
|------|-----------------------|
| Dark | The global dark frame |
| Flat | The global flat field |

---

`> show`

This can be called to view the list of loaded files, and the name under which they are stored.

---
`> move`
`> d`


`> save`

This will prompt a group name, the contents will be saved into `./out/filename.fits`

---

`> terminal`

This [in theory] allows the user to access basic bash commands to navigate/manipulate files/execute other programs, while still inside StarBug. Useful to create input files, move outputs and check results, without losing loaded StarBug files. To exit terminal mode, type `EXIT`, this will bring the user back into normal operation of StarBug.

---

`> clean`

This will delete the contents of `./out/` . To clean all outputs from current directory. **CAREFUL!**

---

`> help`

Prints basic help page.

---

`> exit`

Exits the program.

## Image Reduction

`> build_dark`

Prompted with loaded group name containing dark frames to be combined. This regime takes the mean of all the pixel arrays (and currently just changes the initial file in the group, but in future will create an entirely new fits object). The regime will save the output to the *Dark* preloaded group.

---

`> build_flat`

Prompts loaded group name for group containing flat field frames to be combined. It will use the *Dark* preloaded group as the dark frame to subtract from each frame. If no file exists, then it will not conduct the dark frame subtraction, but will continue to comebine the frames. It scales each frame to a local median, and combines them with pixel-pixel median values. The output is saved into *Flat* preloaded group (again, currently just changes the first instance in list, but finally will create a new file).

---

`> subtract_dark`

Prompts loaded group name, and will iterate through group, subtracting the *Dark* preloaded group from each. Will not conduct regime if no *Dark* is loaded.

---

`> flat_fielding`

Prompts loaded group name, will iterate through group and divide the *Flat* preload from each. If no *Flat* loaded, the, the regime will not be conducted.

`> align`

Calculates the offset between fits images in loaded group

`> exportoffset`

//DEBUGGING PURPOSES
Prompts for filename (if nothing supplied, subtracts .fits off first fitsimage of loaded group as filename), and saves their offsets into file

`> stack`

Prompts offset filename (if nothing supplied, uses offsets stored in loaded group). Function stacks together the images based on their offsets, and crops the resulting data around non-fully fullfilled areas.

## Analysis

`> stats`

Prompts user for loaded group, sigma clip value, iterations. Prints basic stats for each image in group.

## File Manipulation

`> dtype`

Prompts loaded group name, new data type. Converts the fits pixel array to that data type. Currently accepts: *float32 float64* . 

`> header`

Displays header files of each fits image in loaded group

`> update_header`

Prompts loaded group to work on, and header element to alter. Then reads in new value. Loaded group must be saved for change to take effect.
