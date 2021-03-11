all:
	for i in 16; do for lat in scl bcc fcc dia; do python pyrochlore.py $$lat $$i; done; done
