import MySQLdb
import sys

if len(sys.argv) != 5:
  print("<database> <password> <track> <html_file>")
  print("example: python3 loadDescription.py sacCer3 password ensGenesLocal description.txt")
  sys.exit()

database, password, track, html_file = sys.argv[1:5]
html = open(html_file).read()

# input password and organism name
db = MySQLdb.connect(host="localhost",   
                     user="root",        
                     passwd=password,  
                     db=database)       

cur = db.cursor()

sql = 'UPDATE trackDb SET html = %s where tableName=%s'
try:
   # Execute the SQL command
   cur.execute(sql, (html, track))
   # Commit your changes in the database
   db.commit()
except:
   # Rollback in case there is any error
   print("Something went wrong")
   db.rollback()


cur.execute('SELECT tableName, html FROM trackDb where tableName=%s',(track,))
#print all the first cell of all the rows
for row in cur.fetchall():
    print("Value in database:")
    print(row)

db.close()
