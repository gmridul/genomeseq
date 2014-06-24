var casper = require('casper').create();

casper.start('http://www.ncbi.nlm.nih.gov/nuccore/AC189277.4?report=fasta&format=text');

casper.wait(1000, function() {
    this.echo(this.fetchText('div'));
});

casper.run();
