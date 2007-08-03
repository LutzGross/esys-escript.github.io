#!/usr/bin/env python

# $Id: todoBodyMake.py,v 1.4 2005/01/24 07:11:15 paultcochrane Exp $

import xml.dom.minidom
import re
import time

# the head string for the .part page
headString = '<!-- $'
headString += 'Id'
headString += '$ -->\n'
headString += '  <h1>Todo list</h1>\n'
headString += '  <hr class="top" />\n'

# open the todo file
f = open('../../.todo','r')

# parse the document
doc = xml.dom.minidom.parse(f)

# close the todo file (like a good little boy)
f.close()

# grab all of the <note> elements
notes = doc.getElementsByTagName('note')

# regular expression objects for the whitespace at the start and end
# of strings respectively
r1 = re.compile(r'^\s+')
r2 = re.compile(r'\s+$')

# the html string object for the 'done' items
doneHtml = '''
<h2>Completed todo items</h2>
<table class="todo">
  <tr>
    <th class="description">Task description</th>
    <th class="dateAdded">Date added</th>
    <th class="dateCompleted">Date completed</th>
    <th class="comment">Comment</th>
  </tr>
'''

# the html string for the 'todo' items
todoHtml = '''
<h2>Todo items</h2>
<table class="todo">
  <tr>
    <th class="description">Task description</th>
    <th class="dateAdded">Date added</th>
    <th class="priority">Priority</th>
  </tr>
'''

for note in notes:
    # if the thing is done, store what it was, when it was finished, 
    # and the comment
    if note.hasAttribute('done'):
	# grab the attributes
	startTime = note.getAttribute('time')
	startTime = time.strftime('%Y-%m-%d',time.localtime(float(startTime)))
	doneTime = note.getAttribute('done')
	doneTime = time.strftime('%Y-%m-%d',time.localtime(float(doneTime)))
	priority = note.getAttribute('priority')
	# grab the note text
	noteText = note.firstChild.nodeValue
	noteText = r1.sub('',noteText)
	noteText = r2.sub('',noteText)
	# grab the comment text
	comments = note.getElementsByTagName('comment')
	commentNode = comments[0]
	commentText = commentNode.firstChild.nodeValue
	commentText = r1.sub('',commentText)
	commentText = r2.sub('',commentText)

	# now generate the html for 'done' items
	doneHtml += "<tr class=\"done\">\n"
	doneHtml += "  <td class=\"description\">%s</td>\n" % noteText
	doneHtml += "  <td class=\"dateAdded\">%s</td>\n" % startTime
	doneHtml += "  <td class=\"dateCompleted\">%s</td>\n" % doneTime
	doneHtml += "  <td class=\"comment\">%s</td>\n" % commentText
	doneHtml += "</tr>\n"

    else:
	# grab the attributes
	startTime = note.getAttribute('time')
	startTime = time.strftime('%Y-%m-%d',time.localtime(float(startTime)))
	priority = note.getAttribute('priority')
	# grab the note text
	noteText = note.firstChild.nodeValue
	noteText = r1.sub('',noteText)
	noteText = r2.sub('',noteText)
	
	# now generate the html for the 'todo' items
	todoHtml += "<tr class=\"%s\">\n" % priority
	todoHtml += "  <td class=\"description\">%s</td>\n" % noteText
	todoHtml += "  <td class=\"dateAdded\">%s</td>\n" % startTime
	todoHtml += "  <td class=\"priority\">%s</td>\n" % priority
	todoHtml += "</tr>\n"

# finish off the html elements
doneHtml += "</table>\n"
todoHtml += "</table>\n"

# generate the html string to print
todoBodyPartHtml = headString + todoHtml + doneHtml

# write the html to file
f = open("todo_body.part","w")
f.write(todoBodyPartHtml)
f.close()


