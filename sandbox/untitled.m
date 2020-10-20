b = struct();

b.a = 10;

b.b = 'butts';


b.c = [1,2,3,4,5];
b.d = struct();
b.d.butts = 'butt';
b.d.poops = 'poop';






jsonencode(b)


jsondecode('{"a":10,\r\n"b":"butts",\n"c":[1,2,3,4,5],\n"d":{"butts":"butt","poops":"poop"}}')