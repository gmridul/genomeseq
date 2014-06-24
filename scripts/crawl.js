var casper = require('casper').create();

casper.start(casper.cli.get(0));

casper.wait(1000, function() {
    this.echo(this.fetchText('pre'));
});

casper.run();
