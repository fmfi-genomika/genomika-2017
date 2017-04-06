import fileinput
import sys, urllib as ul

for line in fileinput.input():
	sys.stdout.write(ul.unquote_plus(line))
