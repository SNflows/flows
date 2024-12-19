#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# TODO: Imported as is/was from zion:/usr/users/kasoc/flows/ (zion.phys.au.dk: 10.28.0.245),
#       so will need modifications

import argparse
import logging
import os
import datetime
#import getpass
import smtplib
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from astropy.coordinates import Angle
import sys
if '/usr/users/kasoc/Preprocessing/' not in sys.path:
	sys.path.insert(0, '/usr/users/kasoc/Preprocessing/')
#sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), 'flows')))
from kasocutil import psql_connect # 2024-12-19 erik: no trace of this module in '/usr/users/kasoc' ...
#import flows

if __name__ == '__main__':

	# Parse command line arguments:
	parser = argparse.ArgumentParser(description='Send out candidate e-mails.')
	parser.add_argument('-d', '--debug', help='Print debug messages.', action='store_true')
	parser.add_argument('-q', '--quiet', help='Only report warnings and errors.', action='store_true')
	args = parser.parse_args()

	# Set logging level:
	logging_level = logging.INFO
	if args.quiet:
		logging_level = logging.WARNING
	elif args.debug:
		logging_level = logging.DEBUG

	# Setup logging:
	formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
	console = logging.StreamHandler()
	console.setFormatter(formatter)
	logger = logging.getLogger(__name__)
	if not logger.hasHandlers():
		logger.addHandler(console)
	logger.setLevel(logging_level)

	#passwd = getpass.getpass('Password: ')

	with psql_connect('flows_notifier') as conn:
		cursor = conn.cursor()

		cursor.execute("SELECT * FROM flows.targets WHERE target_status='candidate' AND included_in_email IS NULL ORDER BY target_name;")
		results = cursor.fetchall()

		if results:
			logger.info("%d new candidates found.", len(results))

			now = datetime.datetime.utcnow()

			html = """\
			<html>
			<head>
			<style type="text/css">
			body {
				font-family: "Open Sans", Verdana, Arial, sans-serif;
			}
			table.tablesorter {
				background-color: #CDCDCD;
				margin: 10px 0pt 15px;
				width: 100%;
				text-align: left;
			}
			table.tablesorter thead tr th, table.tablesorter tfoot tr th {
				background-color: #e6EEEE;
				border: 1px solid #FFF;
				padding: 4px;
			}
			table.tablesorter thead tr .header {
				background-image: url(/css/tablesorter/bg.gif);
				background-repeat: no-repeat;
				background-position: center right;
				cursor: pointer;
			}
			table.tablesorter tbody td {
				color: Black;
				padding: 4px;
				background-color: #FFF;
				vertical-align: top;
			}
			p.footer {
				font-size: smaller;
			}
			</style>
			</head>
			<body>"""
			html += "<p>Flows, Aarhus University<br>\n" + now.strftime('%d %B %Y %H:%M') + "</p>";
			html += "<p>Dear Flows members,</p>"

			html += "<p>The following candidates have been automatically added to the Flows candidate list:</p>"
			html += '<table class="tablesorter" style="width:100%;"><thead>'
			html += "<tr>"
			html += '<th style="width:10%;">Candidate</th>'
			html += '<th style="width:18%;">RA</th>'
			html += '<th style="width:18%;">Dec</th>'
			html += '<th style="width:18%;">Redshift</th>'
			html += '<th style="width:18%;">Discovery mag.</th>'
			html += '<th style="width:18%;">Discovery date</th>'
			html += "</tr>"
			html += "</thead><tbody>"
			for row in results:
				print(row)

				html += "<tr>"
				html += '<td><a href="https://flows.phys.au.dk/candidates/{0:d}" target="_blank">{1:s}</a></td>'.format(row['targetid'], row['target_name'])
				html += '<td style="text-align:right">{0}</td>'.format(Angle(row['ra']/15, unit='deg').to_string(sep=':', precision=1))
				html += '<td style="text-align:right">{0}</td>'.format(Angle(row['decl'], unit='deg').to_string(sep=':', alwayssign=True, precision=1))
				html += '<td style="text-align:right">{0:.3f}</td>'.format(row['redshift'])
				if row['discovery_mag'] is None:
					html += '<td>&nbsp;</td>'
				else:
					html += '<td style="text-align:right">{0:.2f}</td>'.format(row['discovery_mag'])
				if row['discovery_date'] is None:
					html += '<td>&nbsp;</td>'
				else:
					html += '<td style="text-align:right">{0}</td>'.format(row['discovery_date'].strftime('%Y-%m-%d %H:%M'))
				html += "</tr>\n"
			html += "</tbody></table>"

			html += "<p>Best regards,<br>The Flows Team</p>"
			#html += "<p class=\"footer\">If you no longer wish to receive these e-mails, go to 'My Account' on the TASOC website to disable e-mail notifications about upcoming events.</p>"
			html += "</body>"
			html += "</html>"

			#recipients = ['rasmush@phys.au.dk']
			recipients = ['flows-work.phys@maillist.au.dk']

			# Send an e-mail that the file is ready for download:
			# Create message container - the correct MIME type is multipart/alternative.
			msg = MIMEMultipart('alternative')
			msg['Subject'] = "New Flows candidates"
			msg['From'] = "Flows <rasmush@phys.au.dk>"
			msg['To'] = ', '.join(recipients)

			#msg.attach(MIMEText(text, 'plain'))
			msg.attach(MIMEText(html, 'html'))

			# Send the message via local SMTP server.
			s = smtplib.SMTP('localhost')
			# sendmail function takes 3 arguments: sender's address, recipient's address
			# and message to send - here it is sent as one string.
			s.sendmail(msg['From'], recipients, msg.as_string())
			s.quit()
			logger.info("E-mail sent!")

			for row in results:
				cursor.execute("UPDATE flows.targets SET included_in_email=%s WHERE targetid=%s;", [now, row['targetid']])
			conn.commit()
			logger.info("Targets updated.")

		else:
			logger.warning("No new candidates to broadcast")

		cursor.close()

	logger.info("Done.")
