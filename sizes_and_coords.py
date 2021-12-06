import numpy as np

top19_names = ['Phoenix, AZ', 'Seattle, WA', 'Newark, NJ', 'Las Vegas, NV',
 'Boston, MA', 'Baltimore, MD', 'Los Angeles, CA', 'Minneapolis, MN',
 'Philadelphia, PA', 'San Francisco, CA', 'Houston, TX', 'Chicago, IL', 'Detroit, MI',
 'Charlotte, NC', 'Orlando, FL', 'Atlanta, GA', 'Denver, CO',
 'Dallas/Fort Worth, TX', 'St. Louis, MO']

top19_urb_pops = np.array([2907178,2712338,17800238,1314296,4033503,2075957,11789880,
	2388877,5149499,2995291,3823588,8307561,3903682,760160,1158422,3501725,
	1984585,4146597,2078247])
top19_urb_pops.shape = (19,1)

northern_coords = np.array([
	[53.3856472849378, -1.4816110435792391],[53.52161528235629, -1.1267172945758726],[53.430161199257455, -1.3538710504020839],
	[53.8056666395197, -1.5616251343419858],[53.79113728061139, -1.7527469777986262],[53.64668837224571, -1.7981945778244772],
	[53.47436112230958, -2.2650011085393063],[53.576074929398125, -2.4305221487239774],[53.42567888831289, -2.323574178182925],
	[53.606334307224365, -2.1348334914616998],[53.410810643428306, -2.164506445467592],[53.48734972972736, -2.290170937906172],
	[53.54182966245217, -2.6383251763506124]])

northern_pops = np.array([
	518090,109805,109691,
	474632,349561,162949,
	510746,194189,134022,107926,105878,103886,103608])
northern_pops.shape = (13,1)

northern_names = [
	"Sheffield","Doncaster","Rotherham",
	"Leeds","Bradford","Huddersfield",
	"Manchester","Bolton","Sale","Rochdale","Stockport","Salford","Wigan"]

northern_district_names = ["Bolton","Manchester","Rochdale","Salford","Stockport",
"Trafford","Wigan","Doncaster","Rotherham","Sheffield","Bradford","Kirklees","Leeds"]

northern_district_coords = np.array([
	[34.614     , 85.692     ],[28.76666667, 74.71666667],[36.582     , 69.366     ],
	[29.25      , 77.406     ],[24.636     , 69.45      ],[25.536     , 79.422     ],
	[32.706     , 97.95      ],[31.368     ,  7.71      ],[25.956     , 21.81      ],
	[23.        , 28.        ],[48.        , 45.        ],[38.76      , 47.328     ],
	[47.98333333, 32.95      ]])

northern_district_pops = np.array([194189,510746,107926,
	103886,105878,134022,
	103608,109805,109691,518090,349561,162949,474632])
northern_district_pops.shape = (13,1)

UK_districts = ["Hartlepool","Middlesbrough","Redcar and Cleveland","Stockton-on-Tees",
"Darlington","County Durham","Northumberland","Gateshead","Newcastle upon Tyne",
"North Tyneside","South Tyneside","Sunderland","Halton","Warrington",
"Blackburn with Darwen","Blackpool","Cheshire East","Cheshire West and Chester",
"Allerdale","Barrow-in-Furness","Carlisle","Copeland","Eden","South Lakeland","Burnley",
"Chorley","Fylde","Hyndburn","Lancaster","Pendle","Preston","Ribble Valley",
"Rossendale",
"South Ribble","West Lancashire","Wyre","Bolton","Bury","Manchester","Oldham",
"Rochdale",
"Salford","Stockport","Tameside","Trafford","Wigan","Knowsley","Liverpool","St. Helens",
"Sefton","Wirral","Kingston upon Hull, City of","East Riding of Yorkshire",
"North East Lincolnshire","North Lincolnshire","York","Craven","Hambleton","Harrogate",
"Richmondshire","Ryedale","Scarborough","Selby","Barnsley","Doncaster","Rotherham",
"Sheffield","Bradford","Calderdale","Kirklees","Leeds","Wakefield","Derby","Leicester",
"Rutland","Nottingham","Amber Valley","Bolsover","Chesterfield","Derbyshire Dales",
"Erewash","High Peak","North East Derbyshire","South Derbyshire","Blaby","Charnwood",
"Harborough","Hinckley and Bosworth","Melton","North West Leicestershire",
"Oadby and Wigston","Boston","East Lindsey","Lincoln","North Kesteven","South Holland",
"South Kesteven","West Lindsey","Corby","Daventry","East Northamptonshire","Kettering",
"Northampton","South Northamptonshire","Wellingborough","Ashfield","Bassetlaw",
"Broxtowe",
"Gedling","Mansfield","Newark and Sherwood","Rushcliffe","Herefordshire, County of",
"Telford and Wrekin","Stoke-on-Trent","Shropshire","Cannock Chase","East Staffordshire",
"Lichfield","Newcastle-under-Lyme","South Staffordshire","Stafford",
"Staffordshire Moorlands","Tamworth","North Warwickshire","Nuneaton and Bedworth",
"Rugby","Stratford-on-Avon","Warwick","Bromsgrove","Malvern Hills","Redditch",
"Worcester","Wychavon","Wyre Forest","Birmingham","Coventry","Dudley","Sandwell",
"Solihull","Walsall","Wolverhampton","Peterborough","Luton","Southend-on-Sea",
"Thurrock",
"Bedford","Central Bedfordshire","Cambridge","East Cambridgeshire","Fenland",
"Huntingdonshire","South Cambridgeshire","Basildon","Braintree","Brentwood",
"Castle Point","Chelmsford","Colchester","Epping Forest","Harlow","Maldon","Rochford",
"Tendring","Uttlesford","Broxbourne","Dacorum","East Hertfordshire","Hertsmere",
"North Hertfordshire","St Albans","Stevenage","Three Rivers","Watford",
"Welwyn Hatfield",
"Breckland","Broadland","Great Yarmouth","King's Lynn and West Norfolk",
"North Norfolk",
"Norwich","South Norfolk","Babergh","Forest Heath","Ipswich","Mid Suffolk",
"St Edmundsbury","Suffolk Coastal","Waveney","Barking and Dagenham","Barnet","Bexley",
"Brent","Bromley","Camden","Croydon","Ealing","Enfield","Greenwich","Hackney",
"Hammersmith and Fulham","Haringey","Harrow","Havering","Hillingdon","Hounslow",
"Islington","Kensington and Chelsea","Kingston upon Thames","Lambeth","Lewisham",
"Merton","Newham","Redbridge","Richmond upon Thames","Southwark","Sutton",
"Tower Hamlets","Waltham Forest","Wandsworth","Westminster,City of London","Medway",
"Bracknell Forest","West Berkshire","Reading","Slough","Windsor and Maidenhead",
"Wokingham","Milton Keynes","Brighton and Hove","Portsmouth","Southampton",
"Isle of Wight","Aylesbury Vale","Chiltern","South Bucks","Wycombe","Eastbourne",
"Hastings","Lewes","Rother","Wealden","Basingstoke and Deane","East Hampshire",
"Eastleigh","Fareham","Gosport","Hart","Havant","New Forest","Rushmoor","Test Valley",
"Winchester","Ashford","Canterbury","Dartford","Dover","Gravesham","Maidstone",
"Sevenoaks","Shepway","Swale","Thanet","Tonbridge and Malling","Tunbridge Wells",
"Cherwell","Oxford","South Oxfordshire","Vale of White Horse","West Oxfordshire",
"Elmbridge","Epsom and Ewell","Guildford","Mole Valley","Reigate and Banstead",
"Runnymede","Spelthorne","Surrey Heath","Tandridge","Waverley","Woking","Adur","Arun",
"Chichester","Crawley","Horsham","Mid Sussex","Worthing","Bath and North East Somerset",
"Bristol, City of","Cornwall,Isles of Scilly","Wiltshire","North Somerset",
"South Gloucestershire","Plymouth","Torbay","Bournemouth","Poole","Swindon",
"East Devon",
"Exeter","Mid Devon","North Devon","South Hams","Teignbridge","Torridge","West Devon",
"Christchurch","East Dorset","North Dorset","Purbeck","West Dorset",
"Weymouth and Portland","Cheltenham","Cotswold","Forest of Dean","Gloucester",
"Stroud","Tewkesbury","Mendip","Sedgemoor","South Somerset","Taunton Deane",
"West Somerset","Isle of Anglesey","Gwynedd","Conwy","Denbighshire","Flintshire",
"Wrexham","Ceredigion","Pembrokeshire","Carmarthenshire","Swansea","Neath Port Talbot",
"Bridgend","The Vale of Glamorgan","Cardiff","Rhondda Cynon Taf","Caerphilly",
"Blaenau Gwent","Torfaen","Monmouthshire","Newport","Powys","Merthyr Tydfil",
"Clackmannanshire","Dumfries and Galloway","East Ayrshire","East Lothian",
"East Renfrewshire","Eilean Siar","Falkirk","Fife","Highland","Inverclyde",
"Midlothian","Moray","North Ayrshire","Orkney Islands","Perth and Kinross",
"Scottish Borders","Shetland Islands","South Ayrshire","South Lanarkshire","Stirling",
"Aberdeen City","Aberdeenshire","Argyll and Bute","City of Edinburgh","Renfrewshire",
"West Dunbartonshire","West Lothian","Angus","Dundee City","North Lanarkshire",
"East Dunbartonshire","Glasgow City","Magherafelt","North Down","Dungannon","Belfast",
"Ballymena","Larne","Strabane","Newry and Mourne","Omagh","Armagh","Castlereagh",
"Down","Antrim","Limavady","Derry","Craigavon","Banbridge","Ards","Lisburn",
"Carrickfergus","Ballymoney","Fermanagh","Moyle","Newtownabbey","Cookstown","Coleraine"]

districts_and_pops_dict = {"United Kingdom":63285145,"Great Britain":61470827,
"England and Wales":56170927,"England":53107169,"Wales":3063758,"Scotland":5299900,
"Northern Ireland":1814318,"North East":2596441,"North West":7055961,
"Yorkshire and The Humber":5288212,"East Midlands":4537448,"West Midlands":5608667,
"East":5862418,"London":8204407,"South East":8652784,"South West":5300831,
"Hartlepool":92088,"Middlesbrough":138368,"Redcar and Cleveland":135164,
"Stockton-on-Tees":191824,"Darlington":105584,"Halton":125722,"Warrington":202709,
"Blackburn with Darwen":147657,"Blackpool":142080,"City of Kinston upon Hull":256123,
"East Riding of Yorkshire":334673,"North East Lincolnshire":159735,
"North Lincolnshire":167516,"York":197783,"Derby":248943,"Leicester":329627,
"Rutland":37581,"Nottingham":303899,"Herefordshire":183619,"Telford and Wrekin":166831,
"Stoke-on-Trent":248719,"Bath and North East Somerset":175538,"City of Bristol":428074,
"North Somerset":203091,"South Gloucestershire":263417,"Plymouth":256589,
"Torbay":131193,"Bournemouth":183450,"Poole":148075,"Swindon":209709,
"Peterborough":184457,"Luton":203641,"Southend-on-Sea":174274,"Thurrock":158268,
"Medway":264885,"Bracknell Forest":113696,"West Berkshire":154148,"Reading":155339,
"Slough":140713,"Windsor and Maidenhead":145098,"Wokingham":154943,
"Milton Keynes":249895,"Brighton and Hove":272952,"Portsmouth":205433,
"Southampton":235870,"Isle of Wight":138392,"County Durham":512994,
"Northumberland":316278,"Cheshire East":370736,"Cheshire West and Chester":329526,
"Shropshire":307108,"Cornwall":533760,"Isles of Scilly":2224,"Wiltshire":474319,
"Bedford":157840,"Central Bedfordshire":255644,
"Isle of Anglesey":69913,"Gwynedd":121523,"Conwy":115326,"Denbighshire":93919,
"Flintshire":152666,"Wrexham":135070,"Powys":133071,"Ceredigion":75293,
"Pembrokeshire":122613,"Carmarthenshire":183961,"Swansea":238691,
"Neath Port Talbot":139880,"Bridgend":139410,"The Vale of Glamorgan":126679,
"Cardiff":345442,"Rhondda Cynon Taf":234373,"Merthyr Tydfil":58851,"Caerphilly":178782,
"Blaenau Gwent":69812,"Torfaen":91190,"Monmouthshire":91508,"Newport":145785,
"Greater London":8204407,"Greater Manchester":2685386,"Merseyside":1380770,
"South Yorkshire":1343805,"Tyne And Wear":1104141,"West Midlands":2739733,
"West Yorkshire":2227371,"Buckinghamshire":506550,"Cambridgeshire":622312,
"Cumbria":499817,"Derbyshire":770688,"Devon":747709,"Dorset":413813,
"East Sussex":527209,"Essex":1396599,"Gloucestershire":598289,"Hampshire":1322118,
"Hertfordshire":1119824,"Kent":1466466,"Lancashire":1171558,"Leicestershire":651179,
"Lincolnshire":714768,"Norfolk":859426,"Northamptonshire":693967,
"North Yorkshire":601206,"Nottinghamshire":786796,"Oxfordshire":654791,
"Somerset":531581,"Staffordshire":849546,"Suffolk":730133,"Surrey":1135367,
"Warwickshire":546554,"West Sussex":808919,"Worcestershire":566557,
"City of London":7412,"Barking and Dagenham":187029,"Barnet":357538,"Bexley":232774,
"Brent":312245,"Bromley":310554,"Camden":220087,"Croydon":364815,"Ealing":339314,
"Enfield":313935,"Greenwich":255483,"Hackney":247182,"Hammersmith and Fulham":182445,
"Haringey":255540,"Harrow":240499,"Havering":237927,"Hillingdon":275499,
"Hounslow":254927,"Islington":206285,"Kensington and Chelsea":158251,
"Kingston upon Thames":160436,"Lambeth":304481,"Lewisham":276938,"Merton":200543,
"Newham":310460,"Redbridge":281395,"Richmond upon Thames":187527,"Southwark":288717,
"Sutton":191123,"Tower Hamlets":256012,"Waltham Forest":259742,"Wandsworth":307710,
"Westminster":219582,
"Bolton":277296,"Bury":185422,"Manchester":502902,"Oldham":225157,"Rochdale":211929,
"Salford":234487,"Stockport":283253,"Tameside":219727,"Trafford":227091,"Wigan":318122,
"Knowsley":145903,"Liverpool":465656,"St. Helens":175405,"Sefton":273969,
"Wirral":319837,
"Barnsley":231865,"Doncaster":302468,"Rotherham":257716,"Sheffield":551756,
"Gateshead":200349,"Newcastle upon Tyne":279092,"North Tyneside":201206,
"South Tyneside":148164,"Sunderland":275330,
"Birmingham":1074283,"Coventry":316915,"Dudley":313261,"Sandwell":309042,
"Solihull":206856,"Walsall":269524,"Wolverhampton":249852,"Bradford":523115,
"Calderdale":204170,"Kirklees":422970,"Leeds":750683,"Wakefield":326433,
"Aylesbury Vale":174880,"Chiltern":92652,"South Bucks":67060,"Wycombe":171958,
"Cambridge":122725,"East Cambridgeshire":84245,"Fenland":95461,
"Huntingdonshire":170039,"South Cambridgeshire":149842,
"Allerdale":96444,"Barrow-in-Furness":69056,"Carlisle":107475,"Copeland":70627,
"Eden":52502,"South Lakeland":103713,
"Amber Valley":122521,"Bolsover":76029,"Chesterfield":103788,"Derbyshire Dales":71104,
"Erewash":112249,"High Peak":90982,"North East Derbyshire":99100,
"South Derbyshire":94915,
"East Devon":133272,"Exeter":117063,"Mid Devon":77936,"North Devon":93976,
"South Hams":83563,"Teignbridge":124271,"Torridge":63973,"West Devon":53655,
"Christchurch":47916,"East Dorset":87301,"North Dorset":69002,"Purbeck":45184,
"West Dorset":99275,"Weymouth and Portland":65135,
"Eastbourne":99308,"Hastings":90173,"Lewes":97584,"Rother":90729,"Wealden":149415,
"Basildon":174971,"Braintree":147514,"Brentwood":73841,"Castle Point":87964,
"Chelmsford":168491,"Colchester":173614,"Epping Forest":124880,"Harlow":82177,
"Maldon":61720,"Rochford":83333,"Tendring":138062,"Uttlesford":80032,
"Cheltenham":115645,"Cotswold":83180,"Forest of Dean":82200,"Gloucester":121921,
"Stroud":113074,"Tewkesbury":82269,
"Basingstoke and Deane":168550,"East Hampshire":116010,"Eastleigh":125852,
"Fareham":111931,"Gosport":82669,"Hart":91662,"Havant":120783,"New Forest":176789,
"Rushmoor":94354,"Test Valley":116698,"Winchester":116820,
"Broxbourne":93702,"Dacorum":145298,"East Hertfordshire":138155,"Hertsmere":100379,
"North Hertfordshire":127494,"St Albans":141248,"Stevenage":84247,
"Three Rivers":87921,"Watford":90653,"Welwyn Hatfield":110727,
"Ashford":118405,"Canterbury":150600,"Dartford":97604,"Dover":111718,
"Gravesham":101766,"Maidstone":155764,"Sevenoaks":115351,"Shepway":108199,
"Swale":136324,"Thanet":134402,"Tonbridge and Malling":121087,"Tunbridge Wells":115246,
"Burnley":87032,"Chorley":107591,"Fylde":76098,"Hyndburn":80549,"Lancaster":137823,
"Pendle":89576,"Preston":140054,"Ribble Valley":57292,"Rossendale":68053,
"South Ribble":109181,"West Lancashire":110617,"Wyre":107692,
"Blaby":94132,"Charnwood":165876,"Harborough":85699,"Hinckley and Bosworth":105328,
"Melton":50495,"North West Leicestershire":93670,"Oadby and Wigston":55979,
"Boston":64615,"East Lindsey":136683,"Lincoln":93085,"North Kesteven":108518,
"South Holland":88390,"South Kesteven":134125,"West Lindsey":89352,
"Breckland":131009,"Broadland":124740,"Great Yarmouth":97424,
"King's Lynn and West Norfolk":147936,"North Norfolk":101664,"Norwich":132158,
"South Norfolk":124495,"Corby":61607,"Daventry":78070,"East Northamptonshire":86869,
"Kettering":93846,"Northampton":212492,"South Northamptonshire":85446,
"Wellingborough":75637,"Craven":55459,"Hambleton":89602,"Harrogate":158683,
"Richmondshire":53287,"Ryedale":51893,"Scarborough":108735,"Selby":83547,
"Ashfield":119522,"Bassetlaw":113003,"Broxtowe":109749,"Gedling":113741,
"Mansfield":104551,"Newark and Sherwood":114982,"Rushcliffe":111248,
"Cherwell":142252,"Oxford":150245,"South Oxfordshire":134961,
"Vale of White Horse":121891,"West Oxfordshire":105442,
"Mendip":109406,"Sedgemoor":114919,"South Somerset":162113,"Taunton Deane":110555,
"West Somerset":34588,"Cannock Chase":97582,"East Staffordshire":113858,
"Lichfield":100911,"Newcastle-under-Lyme":123878,"South Staffordshire":108318,
"Stafford":130895,"Staffordshire Moorlands":97209,"Tamworth":76895,
"Babergh":87901,"Forest Heath":60038,"Ipswich":133729,"Mid Suffolk":97076,
"St Edmundsbury":111443,"Suffolk Coastal":124590,"Waveney":115356,
"Elmbridge":131428,"Epsom and Ewell":75191,"Guildford":137580,"Mole Valley":85637,
"Reigate and Banstead":138375,"Runnymede":80501,"Spelthorne":95852,"Surrey Heath":86378,
"Tandridge":83178,"Waverley":121754,"Woking":99493,
"North Warwickshire":62089,"Nuneaton and Bedworth":125409,"Rugby":100496,
"Stratford-on-Avon":120824,"Warwick":137736,
"Adur":61334,"Arun":149811,"Chichester":113995,"Crawley":107053,"Horsham":131540,
"Mid Sussex":140188,"Worthing":104998,"Bromsgrove":93732,"Malvern Hills":74706,
"Redditch":84318,"Worcester":98679,"Wychavon":117074,"Wyre Forest":98048,
"Aberdeen City":222460,"Aberdeenshire":253650,"Angus":116200,"Argyll and Bute":88930,
"Clackmannanshire":51500,"Dumfries and Galloway":151410,"Dundee City":147200,
"East Ayrshire":122690,"East Dunbartonshire":105000,"East Lothian":99920,
"East Renfrewshire":90810,"City of Edinburgh":477940,"Eilean Siar":27690,
"Falkirk":156250,"Fife":365300,"Glasgow City":593060,"Highland":232730,
"Inverclyde":81220,"Midlothian":83450,"Moray":93470,"North Ayrshire":138090,
"North Lanarkshire":337720,"Orkney Islands":21420,"Perth and Kinross":146850,
"Renfrewshire":174700,"Scottish Borders":113880,"Shetland Islands":23240,
"South Ayrshire":112980,"South Lanarkshire":313900,"Stirling":90330,
"West Dunbartonshire":90610,"West Lothian":175300,
"Antrim":53632,"Ards":78047,"Armagh":59651,"Ballymena":64127,"Ballymoney":31276,
"Banbridge":48333,"Belfast":280922,"Carrickfergus":39102,"Castlereagh":67418,
"Coleraine":58959,"Cookstown":37098,"Craigavon":93334,"Derry":108261,"Down":69934,
"Dungannon":58100,"Fermanagh":62006,"Larne":32136,"Limavady":33610,"Lisburn":120486,
"Magherafelt":45132,"Moyle":17062,"Newry and Mourne":100003,"Newtownabbey":85019,
"North Down":79245,"Omagh":51495,"Strabane":39930}