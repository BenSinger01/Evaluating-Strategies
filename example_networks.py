import numpy as np
import sizes_and_coords as sc
import csv

UK_districts = ["Hartlepool","Middlesbrough","Redcar and Cleveland","Stockton-on-Tees","Darlington","County Durham","Northumberland","Gateshead","Newcastle upon Tyne","North Tyneside","South Tyneside","Sunderland","Halton","Warrington","Blackburn with Darwen","Blackpool","Cheshire East","Cheshire West and Chester","Allerdale","Barrow-in-Furness","Carlisle","Copeland","Eden","South Lakeland","Burnley","Chorley","Fylde","Hyndburn","Lancaster","Pendle","Preston","Ribble Valley","Rossendale","South Ribble","West Lancashire","Wyre","Bolton","Bury","Manchester","Oldham","Rochdale","Salford","Stockport","Tameside","Trafford","Wigan","Knowsley","Liverpool","St. Helens","Sefton","Wirral","Kingston upon Hull, City of","East Riding of Yorkshire","North East Lincolnshire","North Lincolnshire","York","Craven","Hambleton","Harrogate","Richmondshire","Ryedale","Scarborough","Selby","Barnsley","Doncaster","Rotherham","Sheffield","Bradford","Calderdale","Kirklees","Leeds","Wakefield","Derby","Leicester","Rutland","Nottingham","Amber Valley","Bolsover","Chesterfield","Derbyshire Dales","Erewash","High Peak","North East Derbyshire","South Derbyshire","Blaby","Charnwood","Harborough","Hinckley and Bosworth","Melton","North West Leicestershire","Oadby and Wigston","Boston","East Lindsey","Lincoln","North Kesteven","South Holland","South Kesteven","West Lindsey","Corby","Daventry","East Northamptonshire","Kettering","Northampton","South Northamptonshire","Wellingborough","Ashfield","Bassetlaw","Broxtowe","Gedling","Mansfield","Newark and Sherwood","Rushcliffe","Herefordshire, County of","Telford and Wrekin","Stoke-on-Trent","Shropshire","Cannock Chase","East Staffordshire","Lichfield","Newcastle-under-Lyme","South Staffordshire","Stafford","Staffordshire Moorlands","Tamworth","North Warwickshire","Nuneaton and Bedworth","Rugby","Stratford-on-Avon","Warwick","Bromsgrove","Malvern Hills","Redditch","Worcester","Wychavon","Wyre Forest","Birmingham","Coventry","Dudley","Sandwell","Solihull","Walsall","Wolverhampton","Peterborough","Luton","Southend-on-Sea","Thurrock","Bedford","Central Bedfordshire","Cambridge","East Cambridgeshire","Fenland","Huntingdonshire","South Cambridgeshire","Basildon","Braintree","Brentwood","Castle Point","Chelmsford","Colchester","Epping Forest","Harlow","Maldon","Rochford","Tendring","Uttlesford","Broxbourne","Dacorum","East Hertfordshire","Hertsmere","North Hertfordshire","St Albans","Stevenage","Three Rivers","Watford","Welwyn Hatfield","Breckland","Broadland","Great Yarmouth","King's Lynn and West Norfolk","North Norfolk","Norwich","South Norfolk","Babergh","Forest Heath","Ipswich","Mid Suffolk","St Edmundsbury","Suffolk Coastal","Waveney","Barking and Dagenham","Barnet","Bexley","Brent","Bromley","Camden","Croydon","Ealing","Enfield","Greenwich","Hackney","Hammersmith and Fulham","Haringey","Harrow","Havering","Hillingdon","Hounslow","Islington","Kensington and Chelsea","Kingston upon Thames","Lambeth","Lewisham","Merton","Newham","Redbridge","Richmond upon Thames","Southwark","Sutton","Tower Hamlets","Waltham Forest","Wandsworth","Westminster,City of London","Medway","Bracknell Forest","West Berkshire","Reading","Slough","Windsor and Maidenhead","Wokingham","Milton Keynes","Brighton and Hove","Portsmouth","Southampton","Isle of Wight","Aylesbury Vale","Chiltern","South Bucks","Wycombe","Eastbourne","Hastings","Lewes","Rother","Wealden","Basingstoke and Deane","East Hampshire","Eastleigh","Fareham","Gosport","Hart","Havant","New Forest","Rushmoor","Test Valley","Winchester","Ashford","Canterbury","Dartford","Dover","Gravesham","Maidstone","Sevenoaks","Shepway","Swale","Thanet","Tonbridge and Malling","Tunbridge Wells","Cherwell","Oxford","South Oxfordshire","Vale of White Horse","West Oxfordshire","Elmbridge","Epsom and Ewell","Guildford","Mole Valley","Reigate and Banstead","Runnymede","Spelthorne","Surrey Heath","Tandridge","Waverley","Woking","Adur","Arun","Chichester","Crawley","Horsham","Mid Sussex","Worthing","Bath and North East Somerset","Bristol, City of","Cornwall,Isles of Scilly","Wiltshire","North Somerset","South Gloucestershire","Plymouth","Torbay","Bournemouth","Poole","Swindon","East Devon","Exeter","Mid Devon","North Devon","South Hams","Teignbridge","Torridge","West Devon","Christchurch","East Dorset","North Dorset","Purbeck","West Dorset","Weymouth and Portland","Cheltenham","Cotswold","Forest of Dean","Gloucester","Stroud","Tewkesbury","Mendip","Sedgemoor","South Somerset","Taunton Deane","West Somerset","Isle of Anglesey","Gwynedd","Conwy","Denbighshire","Flintshire","Wrexham","Ceredigion","Pembrokeshire","Carmarthenshire","Swansea","Neath Port Talbot","Bridgend","The Vale of Glamorgan","Cardiff","Rhondda Cynon Taf","Caerphilly","Blaenau Gwent","Torfaen","Monmouthshire","Newport","Powys","Merthyr Tydfil","Clackmannanshire","Dumfries and Galloway","East Ayrshire","East Lothian","East Renfrewshire","Eilean Siar","Falkirk","Fife","Highland","Inverclyde","Midlothian","Moray","North Ayrshire","Orkney Islands","Perth and Kinross","Scottish Borders","Shetland Islands","South Ayrshire","South Lanarkshire","Stirling","Aberdeen City","Aberdeenshire","Argyll and Bute","City of Edinburgh","Renfrewshire","West Dunbartonshire","West Lothian","Angus","Dundee City","North Lanarkshire","East Dunbartonshire","Glasgow City","Magherafelt","North Down","Dungannon","Belfast","Ballymena","Larne","Strabane","Newry and Mourne","Omagh","Armagh","Castlereagh","Down","Antrim","Limavady","Derry","Craigavon","Banbridge","Ards","Lisburn","Carrickfergus","Ballymoney","Fermanagh","Moyle","Newtownabbey","Cookstown","Coleraine"]

def network(name,scale,power=2,pops=None,coords=None,p=0.5,seed=1,branches=2,width=None,
	join_strength=1):	
	if name=="commuter districts":
		n=404
		A = np.loadtxt("Census_Commute_Data_2011_by_District.csv",delimiter=";")
		A *= 1-np.eye(n)
		A *= scale
	elif name=="commuter northern rate":
		_,A_full = network("commuter districts",scale)
		n=13
		names = ['Sheffield', 'Doncaster', 'Rotherham', 'Leeds', 'Bradford',
		 'Kirklees', 'Manchester', 'Bolton', 'Trafford', 'Rochdale', 'Stockport',
		  'Salford', 'Wigan']
		town_indices = [np.where(np.array(UK_districts)==name)[0] for name in names]
		A_reduced = np.array([[A_full[i,j] for j in town_indices] for i in town_indices])
		A_reduced.shape = (len(names),len(names))
		town_pops = sc.northern_pops
		A = A_reduced/town_pops
	elif name=='air us top 19 urb rate':
		n = 19
		A_two_Chi = np.loadtxt('air_us_top20_traffic_ordered_by_index.csv')
		A_two_Chi[11,:] += A_two_Chi[16,:]
		A_two_Chi[:,11] += A_two_Chi[:,16]
		A_two_Chi = np.delete(A_two_Chi,16,0)
		A = np.delete(A_two_Chi,16,1)
		A = A/sc.top19_urb_pops
		A *= 1 - np.eye(n)
		A *= scale
	elif name=='gravity':
		n = len(pops)
		A = np.zeros((n,n))
		for loc1 in range(n):
			for loc2 in np.delete(range(n),loc1):
				A[loc1,loc2] = pops[loc2]/\
				(np.linalg.norm(coords[loc1]-coords[loc2])**power)
		A *= scale
	return(n,A)

def joined_network(name,scale,power=2,pops=None,coords=None,p=0.5,seed=1,branches=2,
	width=None,join_indices=[[0,17],],join_strength=1):
	n = len(pops)
	A = np.zeros((n,n))
	_,A_part1 = network(name,scale,power=power,coords=coords,p=p,seed=seed,
		branches=branches,width=width,join_strength=join_strength,
		pops=pops[0:n//2])
	_,A_part2 = network(name,scale,power=power,coords=coords,p=p,seed=seed,
		branches=branches,width=width,join_strength=join_strength,
		pops=pops[n//2:n])
	A[0:n//2,0:n//2] = A_part1
	A[n//2:n,n//2:n] = A_part2
	n_joins = len(join_indices)
	if not hasattr(join_strength,'__iter__'):
		join_strength = np.ones(n_joins)*join_strength
	for i in range(n_joins):
		index_pair = join_indices[i]
		A[index_pair[0],index_pair[1]] = A[index_pair[1],index_pair[0]] =\
		 scale*join_strength[i]
	return(n,A)